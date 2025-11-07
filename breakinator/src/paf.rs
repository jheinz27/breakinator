use crate::cli::Cli;
use std::{cmp::min, fs::File, io::{self, BufRead, BufReader, BufWriter, Write}, 
iter::Peekable, process};
use crate::{Breakpoint, classify_break, read_level_class, print_report, print_table};


pub fn process_paf(args: &Cli) ->  Result<(), Box<dyn std::error::Error>>  {
    // check ranges are appropriate
    if !(0.0..=1.0).contains(&args.margin) {
        eprintln!( "error: `--margin {}` is out of range; must be between 0.0 and 1.0 (inclusive)", args.margin );
        process::exit(1);
    }

    let file = File::open(&args.input)?;
    let mut reader = BufReader::new(file).lines().peekable();

    let output = File::create(&args.out)?;
    let mut writer = BufWriter::new(output);
    if args.rcoord {
        writeln!(writer, "#Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification\tbreak1_read\tbreak2_read")?;
    } else { 
        writeln!(writer, "#Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification")?;
    }
    // instantiate empty cluster
    let mut cluster = Vec::with_capacity(5);

    //track read and break level classifications
    let mut read_counts: Vec<u64> = vec![0; 3]; //[fold, chim, pass, unique] 
    let mut break_counts: Vec<u64> = vec![0; 3]; //[fold, chim, pass] 
    let mut reads_pass_filter: u64 = 0; 
    
    while let Some(Ok(_)) = reader.peek() {
        //get all primary and supplementary alignments of a read
        get_clusters(&mut reader, &mut cluster)?;
        //filter read seqments

        let mut filtered = filter_alignments(&cluster, args.min_mapq, args.min_map_len); 
        
        let num_pass = filtered.len(); 
        if num_pass > 0 {
            reads_pass_filter += 1; 
        } 
        
        if num_pass > 1 {
            //get all breaks in a read if there is more than one alignment 
            let all_breaks = determine_break(&mut filtered, &args, &mut read_counts, &mut break_counts); 
            for b in all_breaks {
                writeln!(writer, "{}", b.as_tsv(&args.rcoord))?; 
            }
            } 
        }
    //write results as tsv or print summary  to terminal 
    if args.tabular {
        print_table(reads_pass_filter, &read_counts, &break_counts, args.input.to_string()).expect("error writing to stdout")
    }else {
        print_report(reads_pass_filter, &read_counts, &break_counts,  &args).expect("error writing to stdout"); 
    }
    Ok(())
}

//function to move ahead one read group at a time and update cluster with the group
fn get_clusters<I>(lines: &mut Peekable<I>, cluster: &mut Vec<String>)-> io::Result<()>
where
    I: Iterator<Item= Result<String, std::io::Error>>,
{
    cluster.clear();
    if let Some(Ok(line)) = lines.next() {
        let cur_id = line.split_once('\t').unwrap().0.to_string();
        cluster.push(line);
        //continue grouping lines with same read ID
        while let Some(Ok(next)) = lines.peek() {
            let next_id = next.split_once('\t').unwrap().0;
            if  next_id == cur_id {
                cluster.push(lines.next().unwrap()?)
            } else {
                break;
            }
        }
    };
    Ok(())
}

//filter out read alignments that fail length or mapQ filters or is secondary alignment 
fn filter_alignments(all_maps: &Vec<String>, mapq: u8, map_len: u32) -> Vec<Vec<&str>>{
    let mut passed_filter: Vec<Vec<&str>> =  Vec::new(); 
    for map in all_maps{
        let fields: Vec<&str> = map.split('\t').collect();
        if fields[11].parse::<u8>().expect("MAPQ field was not a valid integer") >= mapq && (fields[3].parse::<u32>().expect("Alignment end not a valid integer")- fields[2].parse::<u32>().expect("Alignment start not a valid integer")) >= map_len {
            passed_filter.push(fields);
        }
    }
    return passed_filter; 
}

//find the breakpoints and a read that has at least one supplementary alignment
fn determine_break(clust: &mut Vec<Vec<&str>>, args: &Cli, read_counts: &mut Vec<u64>, break_counts: &mut Vec<u64>) ->  Vec<Breakpoint> {
    //sort by start location of aligment in read 
    clust.sort_by_key(|line| {line[2].parse::<i32>().expect("invalid int in field 3")});
    let mut out: Vec<Breakpoint> = Vec::new(); 
    let mut labels: Vec<u32> = vec![0; 3];  //[fold, chim, pass count]

    //check every concurrent alignment in a read with n split alignments 
    for i in 0..(clust.len() - 1) {
        //get strands of both alignments 
        let mut directions = String::new(); 
   
        //first break
        let b1 = match clust[i][4]{
            //First break on forward strand, take end of first alignment
            "+" => {
                directions.push('>'); 
                vec![clust[i][5], clust[i][8]]}
            //First break on reverse strand, take start of first alignment
            "-" => { directions.push('<'); 
                vec![clust[i][5], clust[i][7]]}
            // handle unexpected strand symbol
            _ => {
                eprintln!("Warning: unexpected strand symbol {}", clust[i][4]);
                directions.push('?');
                vec![clust[i][5], clust[i][8]] 
                }         
        };
    
        let b2 = match clust[i+1][4] {
            //second break on forward strand, take start of later alignment 
            "+" => {directions.push('>'); vec![clust[i+1][5], clust[i+1][7]]} 
            //second break on reverse strand, take end of later alignment
            "-" => {directions.push('<'); vec![clust[i+1][5], clust[i+1][8]]}
            // handle unexpected strand symbol
            _ => {
                eprintln!("Warning: unexpected strand symbol {}", clust[i][4]);
                directions.push('?');
                vec![clust[i][5], clust[i][8]] 
                }     
           
        }; 
        
        let read_length = clust[i][1].parse::<u32>().expect("invalid int in field 2");  

        //set mapq val to be min of mapq value for both sides of breakpoint
        let mapq = min(clust[i][11].parse::<u8>().expect("invalid int in field 11"), clust[i+1][11].parse::<u8>().expect("invalid int in field 11"));
        //let mut break_info = vec![b1[0].to_string(), b1[1].to_string(), directions, b2[0].to_string(), b2[1].to_string(),mapq.to_string(), clust[0][0].to_string()];  

        let mut break_info = Breakpoint{b1_chr: b1[0].to_string() ,b1_loc: b1[1].parse::<i64>().expect("Failed to parse b1_loc as i64"), directions: directions,
            b2_chr: b2[0].to_string(), b2_loc: b2[1].parse::<i64>().expect("Failed to parse b2_loc as i64"),
            mapq: mapq, read_id: clust[0][0].to_string(), read_len: read_length, label: None, 
            r1_loc: clust[i][3].parse::<u32>().expect("Failed to parse r1_loc as u32"),
            r2_loc: clust[i+1][2].parse::<u32>().expect("Failed to parse r2_loc as u32")};  
       
        
        //let rlen = clust[i][1].parse::<f32>().expect("invalid int in field 2"); 
        //let rbreak = clust[i][3].parse::<f32>().expect("invalid int in field 4"); 

        //get artifact or pass classification 
        //let label = classify_break( &break_info, args, rlen, rbreak, break_counts); 
        let label = classify_break( &break_info, args, break_counts);  
        match label.as_str() {
            "Foldback" => labels[0] += 1,
            "Chimeric" => labels[1] += 1,
            _ => labels[2] += 1,
        }
        //break_info.push(label); 
        break_info.label = Some(label); 
    
        out.push(break_info); 
    
    }
    //classify the read based on the classification of all breakpoints on the read
    read_level_class(&labels, read_counts); 
    return out; 
    
} 




