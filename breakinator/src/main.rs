use std::{fs::File, io::{BufRead, BufReader, BufWriter, Write},process};
use clap::Parser;
use breakinator::{Cli, get_clusters, filter_alignments, determine_break, print_report, print_table};


fn main() -> Result<(), Box<dyn std::error::Error>> {
    //read in args 
    let args = Cli::parse();
    // check ranges are appropriate
    if !(0.0..=1.0).contains(&args.margin) {
        eprintln!( "error: `--margin {}` is out of range; must be between 0.0 and 1.0 (inclusive)", args.margin );
        process::exit(1);
    }

    let file = File::open(&args.input)?;
    let mut reader = BufReader::new(file).lines().peekable();

    let output = File::create(&args.out)?;
    let mut writer = BufWriter::new(output);
    writeln!(writer, "Break1_chr\tBreak1_loc\tBreak_direction\tBreak2_chr\tBreak2_loc\tMapQ\tRead_ID\tClassification")?; 

    // instantiate empty cluster
    let mut cluster = Vec::with_capacity(4);

    //track read and break level classifications
    let mut read_counts: Vec<i64> = vec![0; 3]; //[fold, chim, pass, unique] 
    let mut break_counts: Vec<i64> = vec![0; 3]; //[fold, chim, pass] 
    let mut reads_pass_filter: i64 = 0; 
    
    while let Some(Ok(_)) = reader.peek() {
        //get all primary and supplementary alignments of a read
        let _ = get_clusters(&mut reader, &mut cluster);
        //filter read seqments
        let mut filtered = filter_alignments(&cluster, args.min_mapq, args.min_map_len); 
        
        let num_pass = filtered.len(); 
        if num_pass > 0 {
            reads_pass_filter += 1; 
        } 
        if num_pass > 1 {
            //get all breaks in a read if there is more than one alignment 
            let breaks = determine_break(&mut filtered, &args, &mut read_counts, &mut break_counts)?; 
            for r in breaks {
                writeln!(writer, "{}", r.join("\t"))?; 
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
