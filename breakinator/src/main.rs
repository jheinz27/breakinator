use clap::Parser;
use breakinator::{Cli, paf, sam};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    //read in args 
    let args = Cli::parse();
    
    if args.paf { 
        if !args.input.to_lowercase().ends_with(".paf") {
            eprintln!("ERROR: File name does not end in paf-ensure file is paf");
        }   
        paf::process_paf(&args)?;
    } else {
        if args.input.to_lowercase().ends_with(".paf") {
            eprintln!("ERROR: include --paf argument if running on paf file");
            std::process::exit(1);
        }   
        
        let is_cram = args.input.to_lowercase().ends_with(".cram");
        if is_cram && args.genome.is_none() {
            eprintln!("ERROR: --genome <FASTA> is required when reading CRAM files.");
            std::process::exit(1);
        }
        sam::process_sam(&args, is_cram)?; 
    }
    Ok(())
}
 