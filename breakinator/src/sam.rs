use rust_htslib::{
    bam::{self, Read, Record,},
    errors::Error as BamError,};
use crate::cli::Cli;
use std::{collections::HashMap, iter::Peekable};
use rust_htslib::{bam::record::{Cigar}};
use rust_htslib::bam::ext::BamRecordExtensions;
use std::cmp::min;
use rust_htslib::bam::HeaderView;
use std::{fs::File, io::{BufWriter, Write}};
use crate::{Breakpoint, classify_break, read_level_class, print_report, print_table};

//alignment info needed to call breakpoints, kept small so many can be held in memory
struct Aln {
    tid: i32,
    pos: i64,
    ref_end: i64,
    is_reverse: bool,
    mapq: u8,
    q_start: u32,
    q_end_clip: u32,
}

impl Aln {
    fn from_record(rec: &Record) -> Aln {
        Aln { tid: rec.tid(), pos: rec.pos(), ref_end: rec.reference_end(), is_reverse: rec.is_reverse(),
            mapq: rec.mapq(), q_start: query_loc(rec), q_end_clip: query_end(rec) }
    }
}

pub fn process_sam(args: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    // read in file, "-" reads from stdin
    let mut sam_reader = if args.input == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(&args.input)?
    };
    sam_reader.set_threads(args.threads)?;
    // attach reference if given (needed for CRAM, ignored otherwise)
    if let Some(genome) = &args.genome {
        sam_reader.set_reference(genome)?;
    }

    let header = sam_reader.header().to_owned();
    // warn if input is position sorted, reads must be grouped by read ID unless running with --sorted
    let is_pos_sorted = String::from_utf8_lossy(header.as_bytes()).lines()
        .any(|l| l.starts_with("@HD") && l.split('\t').any(|f| f == "SO:coordinate"));
    if is_pos_sorted && !args.sorted {
        eprintln!("WARNING: input appears to be sorted by position. Use -s/--sorted, or group reads by read ID (e.g. samtools collate or samtools sort -n)");
    }
    //let header: Header = Header::from_template(&hv);
    //et mut out = Writer::from_stdout(&header, Format::Sam)?;
    let output = File::create(&args.out)?;
    let mut writer = BufWriter::new(output);
    if args.rcoord {
        writeln!(writer, "#Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification\tbreak1_read\tbreak2_read")?;
    } else {
        writeln!(writer, "#Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification")?;
    }

    //track read and break level classifications
    let mut read_counts: Vec<u64> = vec![0; 3]; //[fold, chim, pass, unique]
    let mut break_counts: Vec<u64> = vec![0; 3]; //[fold, chim, pass]

    let reads_pass_filter = if args.sorted {
        process_sorted(&mut sam_reader, args, &header, &mut writer, &mut read_counts, &mut break_counts)?
    } else {
        process_grouped(&mut sam_reader, args, &header, &mut writer, &mut read_counts, &mut break_counts)?
    };

    if args.tabular {
        print_table(reads_pass_filter, &read_counts, &break_counts, args.input.to_string()).expect("error writing to stdout")
    }else {
        print_report(reads_pass_filter, &read_counts, &break_counts,  &args).expect("error writing to stdout");
    }

    Ok(())
}

//input grouped by read ID: process one read group at a time, returns number of reads that passed filter
fn process_grouped(sam_reader: &mut bam::Reader, args: &Cli, header: &HeaderView, writer: &mut impl Write,
    read_counts: &mut Vec<u64>, break_counts: &mut Vec<u64>) -> Result<u64, Box<dyn std::error::Error>> {
    //peakable iterator of file
    let mut sam_iter = sam_reader.records().peekable();

    //vector that store all alignments of one read (cluster of alignments)
    let mut cluster: Vec<Record> = Vec::with_capacity(5);
    let mut reads_pass_filter: u64 = 0;

    while sam_iter.peek().is_some() {
        //move forward by one read group
        if let Err(e) = get_clusters(&mut sam_iter, &mut cluster) {
            read_error_exit(args, e);
        }
        let filtered = filter_alignments(&cluster, &args);

        let num_pass = filtered.len();
        if num_pass > 0 {
            reads_pass_filter += 1;
        }
        if num_pass > 1 {
            let read_len = get_read_len(filtered[0]);
            let read_id = std::str::from_utf8(filtered[0].qname()).unwrap();
            let mut alns: Vec<Aln> = filtered.iter().map(|rec| Aln::from_record(rec)).collect();
            let all_breaks = determine_break(&mut alns, read_id, read_len, &args, read_counts, break_counts, header);
            for b in all_breaks {
                writeln!(writer, "{}", b.as_tsv(&args.rcoord))?;
            }

        }
    }
    Ok(reads_pass_filter)
}

//coordinate sorted input: single pass that keeps only alignments of reads with an SA tag,
//then groups them by read ID in memory. returns number of reads that passed filter
fn process_sorted(sam_reader: &mut bam::Reader, args: &Cli, header: &HeaderView, writer: &mut impl Write,
    read_counts: &mut Vec<u64>, break_counts: &mut Vec<u64>) -> Result<u64, Box<dyn std::error::Error>> {
    //read ID -> (read length, alignments that passed filter)
    let mut split_reads: HashMap<Vec<u8>, (u32, Vec<Aln>)> = HashMap::new();
    let mut reads_pass_filter: u64 = 0;
    let mut num_buffered: u64 = 0;

    let mut rec = Record::new();
    while let Some(result) = sam_reader.read(&mut rec) {
        if let Err(e) = result {
            read_error_exit(args, e);
        }
        if !passes_filter(&rec, args) {
            continue;
        }
        if rec.aux(b"SA").is_err() {
            //no supplementary alignments, read can not contain a breakpoint
            reads_pass_filter += 1;
            continue;
        }
        let entry = split_reads.entry(rec.qname().to_vec()).or_insert_with(|| (get_read_len(&rec), Vec::with_capacity(2)));
        entry.1.push(Aln::from_record(&rec));
        num_buffered += 1;
    }
    eprintln!("Sorted mode: buffered {} alignments from {} reads with SA tags", num_buffered, split_reads.len());

    //sort by read ID so output order is deterministic
    let mut split_reads: Vec<(Vec<u8>, (u32, Vec<Aln>))> = split_reads.into_iter().collect();
    split_reads.sort_unstable_by(|a, b| a.0.cmp(&b.0));

    //SA reads are counted here, once per read, as at least one of their alignments passed filter
    for (qname, (read_len, mut alns)) in split_reads {
        reads_pass_filter += 1;
        if alns.len() > 1 {
            let read_id = std::str::from_utf8(&qname).unwrap();
            let all_breaks = determine_break(&mut alns, read_id, read_len, &args, read_counts, break_counts, header);
            for b in all_breaks {
                writeln!(writer, "{}", b.as_tsv(&args.rcoord))?;
            }
        }
    }
    Ok(reads_pass_filter)
}

//stop on unreadable record instead of silently skipping it
fn read_error_exit(args: &Cli, e: BamError) -> ! {
    eprintln!("ERROR: failed to read alignment record from {}: {}", args.input, e);
    eprintln!("If piping SAM text, make sure the header is included (e.g. samtools view -h)");
    std::process::exit(1);
}

//function to move ahead one read group at a time for SAM/BAM/CRAM
fn get_clusters<I>(records: &mut Peekable<I>, cluster: &mut Vec<Record>)-> Result<(), BamError>
where
    I: Iterator<Item= Result<Record,BamError>>,
{
    cluster.clear();

    if let Some(result) = records.next() {
        //stop on unreadable record instead of silently skipping it
        let record = result?;
        let cur_id = record.qname().to_vec();
        cluster.push(record);

        while let Some(Ok(next)) = records.peek() {
            if next.qname() == cur_id.as_slice() {
                let rec= records.next().unwrap()?;
                cluster.push(rec);
            } else {
                break;
            }
        }
    };
    return Ok(());
}

//filter out read alignments that fail length or mapQ filters or is secondary alignment or unmapped
fn filter_alignments<'a>(all_maps:&'a Vec<Record>, args: &Cli) ->  Vec<&'a Record>{
    return all_maps.iter().filter(|alignment| passes_filter(alignment, args)).collect();
}

//check alignment passes length and mapQ filters and is not secondary alignment or unmapped
fn passes_filter(alignment: &Record, args: &Cli) -> bool {
    if alignment.is_unmapped() || alignment.is_secondary() || alignment.mapq() < args.min_mapq  {
        return false;
    }
    let mut qlen = 0;
    //parse cigar string to determine query length in alignment
    for c in alignment.cigar().iter() {
        match c {
            //consumes query
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l) => {qlen += *l},
            _ => {}
        }
        if qlen >= args.min_map_len {
            return true;
        }
    }
    return false;
}

fn determine_break(clust: &mut Vec<Aln>, read_id: &str, read_length: u32, args: &Cli, read_counts: &mut Vec<u64>, break_counts: &mut Vec<u64>, header:&HeaderView )-> Vec<Breakpoint> {
    //sort by start location of aligment in read
    clust.sort_by_key(|aln| aln.q_start);

    let mut out: Vec<Breakpoint> = Vec::new();
    let mut labels: Vec<u32> = vec![0; 3];  //[fold, chim, pass count]


    //check every concurrent alignment in a read with n split alignments
    for pair in clust.windows(2) {

        let cur = &pair[0];
        let next = &pair[1];

        //get strands of both alignments
        let mut directions = String::new();

        let chr = std::str::from_utf8(header.tid2name(cur.tid as u32)).unwrap();
        let loc;
        //first break
        if cur.is_reverse {
            //First break on reverse strand, take start of first alignment
            directions.push('<');
            loc = cur.pos;
        } else {
            //First break on forward strand, take end of first alignment
            directions.push('>');
            loc = cur.ref_end;
        }

        //second break
        let next_chr = std::str::from_utf8(header.tid2name(next.tid as u32)).unwrap();
        let next_loc;
        if next.is_reverse {
            //second break on reverse strand, take end of later alignment
            directions.push('<');
            next_loc =  next.ref_end;
        } else {
            //second break on forward strand, take start of later alignment
            directions.push('>');
            next_loc = next.pos;
        }

        //set mapq val to be min of mapq value for both sides of breakpoint
        let mapq = min(cur.mapq, next.mapq);

        let mut break_info = Breakpoint{b1_chr: chr.to_string() ,b1_loc: loc, directions: directions,
             b2_chr: next_chr.to_string(), b2_loc: next_loc,
             mapq: mapq, read_id: read_id.to_string(), read_len: read_length, label: None, r1_loc: (read_length - cur.q_end_clip), r2_loc: next.q_start};

        //let rlen = clust[i][1].parse::<f32>().expect("invalid int in field 2");
        //let rbreak = clust[i][3].parse::<f32>().expect("invalid int in field 4");

        //get artifact or pass classification
        let label = classify_break( &break_info, args, break_counts);

        match label.as_str() {
            "Foldback" => labels[0] += 1,
            "Chimeric" => labels[1] += 1,
            _ => labels[2] += 1,
        }

        break_info.label = Some(label);
        out.push(break_info);

    }
    read_level_class(&labels, read_counts);
    return out;
}

fn query_loc(rec: &Record) -> u32 {
    let cigar = rec.cigar();
    let left = cigar.iter().next().and_then(|c| match *c {
        Cigar::SoftClip(l) | Cigar::HardClip(l) => Some(l),
        _ => None,
    }).unwrap_or(0);

    let right = cigar.iter().rev().next().and_then(|c| match *c {
        Cigar::SoftClip(l) | Cigar::HardClip(l) => Some(l),
        _ => None,
    }).unwrap_or(0);

    if rec.is_reverse() {right} else {left}
}

fn query_end(rec: &Record) -> u32 {
    let cigar = rec.cigar();
    let left = cigar.iter().rev().next().and_then(|c| match *c {
        Cigar::SoftClip(l) | Cigar::HardClip(l) => Some(l),
        _ => None,
    }).unwrap_or(0);

    let right = cigar.iter().next().and_then(|c| match *c {
        Cigar::SoftClip(l) | Cigar::HardClip(l) => Some(l),
        _ => None,
    }).unwrap_or(0);

    if rec.is_reverse() {right} else {left}
}


fn get_read_len(rec: &Record) -> u32  {
    let mut qlen = 0;
    //parse cigar string to determine total query length
    for c in rec.cigar().iter() {
        match c {
            //consumes query
            Cigar::SoftClip(l) | Cigar::HardClip(l) | Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l) => {qlen += *l},
            _ => {}
        }
    }
    return qlen;
}
