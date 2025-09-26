use clap::Parser;


#[derive(Parser, Debug)]
#[command( name = "breakinator", about = "Flag foldbacks and chimeric reads from PAF input", version = "1.0")]

pub struct Cli {
    // 
    #[arg(short = 'i', long, value_name = "FILE", required = true, help="PAF file sorted by read IDs")]
    pub input: String,

    // Minimum mapping quality (integer)
    #[arg(short = 'q', long, value_name = "INT", default_value_t = 10, help = "Minimum mapping quality")]
    pub min_mapq: i32,

    // Minimum alignment length (bps)
    #[arg(short = 'a', long,  value_name = "INT", default_value_t = 200, help = "Minimum alignment length (bps)")]
    pub min_map_len: i32,

    // Only report palindromic foldback reads within margin
    #[arg(long, value_name = "BOOL",  default_value_t = false , help="Only report palindromic foldback reads within margin")]
    pub sym: bool,

    // [0-1], With --sym, Proportion from center on either side to be considered foldback artifact
    #[arg(short, long, value_name = "FLOAT", default_value_t = 0.05, help = "[0-1], With --sym, Proportion from center on either side to be considered foldback artifact")]
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
}