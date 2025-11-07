use clap::Parser;


#[derive(Parser, Debug)]
#[command( name = "breakinator", about = "Flag foldbacks and chimeric reads from SAM/BAM/CRAM or PAF input", version = "1.0")]

pub struct Cli {
    // 
    #[arg(short = 'i', long, value_name = "FILE", required = true, help="SAM/BAM/CRAM file sorted by read IDs")]
    pub input: String,

    // input is PAF file
    #[arg(long,value_name = "BOOL", default_value_t = false, help = "Input file is PAF")]
    pub paf: bool,

    // Minimum mapping quality (integer)
    #[arg(short = 'q', long, value_name = "INT", default_value_t = 10, help = "Minimum mapping quality")]
    pub min_mapq: u8,

    // Minimum alignment length (bps)
    #[arg(short = 'a', long,  value_name = "INT", default_value_t = 200, help = "Minimum alignment length (bps)")]
    pub min_map_len: u32,

    // Only report palindromic foldback reads within margin
    #[arg(long, value_name = "BOOL",  default_value_t = false , help="Report all foldback reads, not just those with breakpoint within margin of middle of read")]
    pub no_sym: bool,

    // reference genome used for cram compression
    #[arg(short, long, value_name = "FASTA", help="Reference genome FASTA used (must be provided for CRAM input")]
    pub genome: Option<String>,

    // [0-1], With --no_sym, Proportion from center on either side to be considered foldback artifact
    #[arg(short, long, value_name = "FLOAT", default_value_t = 0.1, help = "[0-1], Proportion from center of read on either side to be considered sym foldback artifact")]
    pub margin: f32,
    
    // Print read coordinates of breakpoint in output
    #[arg(long,value_name = "BOOL", default_value_t = false, help= "Print read coordinates of breakpoint in output" )]
    pub rcoord: bool,

    // Output file name
    #[arg(short = 'o',long, value_name = "FILE", default_value = "breakinator_out.txt", help= "Output file name")]
    pub out: String,

    // Minimum distance to be considered chimeric
    #[arg(short, long, value_name = "INT", default_value_t = 1_000_000, help = "Minimum distance to be considered chimeric")]
    pub chim: i32,

    // Max distance to be considered foldback
    #[arg(short, long, value_name = "INT", default_value_t = 200, help = "Max distance to be considered foldback")]
    pub fold: i32,

    // Return report as a tsv file (useful for evaluating multiple files)
    #[arg(long,value_name = "BOOL", default_value_t = false, help = "Print a TSV table instead of the default report (useful if evaluating multiple samples)")]
    pub tabular: bool,

    // number of threads to use
    #[arg(short, long,value_name = "INT", default_value_t = 2, help = "Number of threads to use for BAM/CRAM I/O")]
    pub threads: usize
}