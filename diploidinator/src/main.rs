use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use diploidinator::*;



fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        eprintln!("ERROR: Too few arguments\nUsage: {} <hap1.paf> <hap2.paf>", args[0]);
        std::process::exit(1);
    }
    let stdout = io::stdout();
    let mut handle = stdout.lock();

    let file1 = File::open(&args[1])?;
    let file2 = File::open(&args[2])?;

    let mut reader1 = BufReader::new(file1).lines().peekable();
    let mut reader2 = BufReader::new(file2).lines().peekable();

    let mut cluster1 = Vec::with_capacity(10); 
    let mut cluster2 = Vec::with_capacity(10); 

    while let Some(Ok(_)) = reader1.peek() {
        let _ = get_clusters(&mut reader1, &mut cluster1);
        let _ = get_clusters(&mut reader2, &mut cluster2);
        if let Some(winner_clust) = compare_clusters(&cluster1, &cluster2) {
            writeln!(handle, "{}", winner_clust.join("\n"))?;
        }
    }
    Ok(())
}