use rust_htslib::{
    bam::{self, Read, Record,},
    errors::Error as BamError,};
use crate::cli::Cli;
use std::{iter::Peekable};
use rust_htslib::{bam::record::{Cigar}};
use rust_htslib::bam::ext::BamRecordExtensions;
use std::cmp::min;
use rust_htslib::bam::HeaderView;
use std::{fs::File, io::{BufWriter, Write}};
use crate::{Breakpoint, classify_break, read_level_class, print_report, print_table};

pub fn process_sam(args: &Cli, is_cram: bool) -> Result<(), Box<dyn std::error::Error>> { 
    // read in file 
    let mut sam_reader = bam::Reader::from_path(&args.input).expect("Failed to open file");
    sam_reader.set_threads(args.threads)?;
    // If CRAM, attach reference
    if is_cram {
        sam_reader.set_reference(args.genome.as_ref().unwrap())?;
    }

    let header = sam_reader.header().to_owned();
    //let header: Header = Header::from_template(&hv); 
    //et mut out = Writer::from_stdout(&header, Format::Sam)?;
    let output = File::create(&args.out)?;
    let mut writer = BufWriter::new(output);
    if args.rcoord {
        writeln!(writer, "#Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification\tbreak1_read\tbreak2_read")?;
    } else { 
        writeln!(writer, "#Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification")?;
    }
 
    //peakable iterator of file 
    let mut sam_iter = sam_reader.records().peekable();

    //vector that store all alignments of one read (cluster of alignments)
    let mut cluster: Vec<Record> = Vec::with_capacity(5); 

    //track read and break level classifications
    let mut read_counts: Vec<u64> = vec![0; 3]; //[fold, chim, pass, unique] 
    let mut break_counts: Vec<u64> = vec![0; 3]; //[fold, chim, pass] 
    let mut reads_pass_filter: u64 = 0; 
    
    while sam_iter.peek().is_some() { 
        //move forward by one read group 
        get_clusters(&mut sam_iter, &mut cluster)?;
        let mut filtered = filter_alignments(&cluster, &args); 
        
        let num_pass = filtered.len(); 
        if num_pass > 0 {
            reads_pass_filter += 1; 
        } 
        if num_pass > 1 {  
            let all_breaks = determine_break(&mut filtered, &args, &mut read_counts, &mut break_counts, &header);
            for b in all_breaks {
                writeln!(writer, "{}", b.as_tsv(&args.rcoord))?; 
            }
            
        }
    }

    if args.tabular {
        print_table(reads_pass_filter, &read_counts, &break_counts, args.input.to_string()).expect("error writing to stdout")
    }else {
        print_report(reads_pass_filter, &read_counts, &break_counts,  &args).expect("error writing to stdout"); 
    }

    Ok(())
}

//function to move ahead one read group at a time for SAM/BAM/CRAM
fn get_clusters<I>(records: &mut Peekable<I>, cluster: &mut Vec<Record>)-> Result<(), BamError>
where 
    I: Iterator<Item= Result<Record,BamError>>,
{
    cluster.clear();

    if let Some(Ok(record)) = records.next() {
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
    let mut passed_filter: Vec<&Record> =  Vec::new(); 
    for alignment in all_maps{
        if alignment.is_unmapped() || alignment.is_secondary() || alignment.mapq() < args.min_mapq  {
            continue; 
        } else {

        }
        let mut is_over_min_len = false ; 
        let mut qlen = 0;
        //parse cigar string to determine query length in alignment
        for c in alignment.cigar().iter() {
            match c {
                //consumes query  
                Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l) => {qlen += *l},
                _ => {}
            }
            if qlen >= args.min_map_len {
                is_over_min_len = true; 
                break; 
            }
        }
        if is_over_min_len {
        
            passed_filter.push(alignment);
        }
    }
    return passed_filter; 
}

fn determine_break(clust: &mut Vec<&Record>, args: &Cli, read_counts: &mut Vec<u64>, break_counts: &mut Vec<u64>, header:&HeaderView )-> Vec<Breakpoint> {
    //sort by start location of aligment in read 
    let read_length = get_read_len(&clust[0]); 

    clust.sort_by_key(|rec| query_loc(rec));
    
    let mut out: Vec<Breakpoint> = Vec::new(); 
    let mut labels: Vec<u32> = vec![0; 3];  //[fold, chim, pass count]
    

    //check every concurrent alignment in a read with n split alignments 
    for pair in clust.windows(2) {

        let cur = &pair[0];
        let next = &pair[1];

        //get strands of both alignments 
        let mut directions = String::new();

        let chr = std::str::from_utf8(header.tid2name(cur.tid() as u32)).unwrap(); 
        let loc; 
        //first break 
        if cur.is_reverse() {
            //First break on reverse strand, take start of first alignment
            directions.push('<'); 
            loc = cur.pos(); 
        } else {
            //First break on forward strand, take end of first alignment
            directions.push('>'); 
            loc = cur.reference_end(); 
        }
        
        //second break 
        let next_chr = std::str::from_utf8(header.tid2name(next.tid() as u32)).unwrap();
        let next_loc; 
        if next.is_reverse() {
            //second break on reverse strand, take end of later alignment 
            directions.push('<'); 
            next_loc =  next.reference_end(); 
        } else {
            //second break on forward strand, take start of later alignment 
            directions.push('>'); 
            next_loc = next.pos();
        }

        //set mapq val to be min of mapq value for both sides of breakpoint
        let mapq = min(cur.mapq(), next.mapq()); 
        let read_id = std::str::from_utf8(cur.qname()).unwrap().to_string();
        
        let mut break_info = Breakpoint{b1_chr: chr.to_string() ,b1_loc: loc, directions: directions,
             b2_chr: next_chr.to_string(), b2_loc: next_loc,
             mapq: mapq, read_id: read_id, read_len: read_length, label: None, r1_loc: (read_length - query_end(cur)), r2_loc: query_loc(next)};  
        
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