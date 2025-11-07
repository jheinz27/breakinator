pub mod cli;
pub use cli::Cli;
pub mod paf;
pub mod sam; 
use std::{ env, io::{self, Write}};

//let mut break_info = vec![b1.0.to_string(), b1.1.to_string(), directions, b2.0.to_string(), b2.1.to_string(),mapq.to_string(), read_id];  

pub struct Breakpoint {
    pub b1_chr: String, 
    pub b1_loc: i64, 
    pub directions: String, 
    pub b2_chr: String, 
    pub b2_loc: i64, 
    pub mapq: u8, 
    pub read_id: String, 
    pub read_len: u32,
    pub label: Option<String>, 
    pub r1_loc: u32, 
    pub r2_loc: u32, 
}
impl Breakpoint {
    pub fn as_tsv(&self, rcoords:&bool) -> String {
        let mut fields = vec![
            self.b1_chr.clone(),
            self.b1_loc.to_string(),
            self.directions.clone(),
            self.b2_chr.clone(),
            self.b2_loc.to_string(),
            self.mapq.to_string(),
            self.read_id.clone(),
            self.label.clone().unwrap_or_else(|| "NA".to_string()),
        ];

        if *rcoords {
            fields.push(self.r1_loc.to_string()); 
            fields.push(self.r2_loc.to_string()); 
        }
        
        fields.join("\t")

    }
}


//classify break as either chimeric, foldback, or pass
pub fn classify_break(brk: &Breakpoint,  args: &Cli, break_counts: &mut Vec<u64>) -> String {
    let dist = (brk.b2_loc - brk.b1_loc).abs() as i32; 
    if brk.b1_chr != brk.b2_chr{
        break_counts[1] +=1; 
        return String::from("Chimeric"); 
    } else if dist >= args.chim {
        break_counts[1] +=1; 
        return String::from("Chimeric"); 
    } else if (brk.directions == "<>" || brk.directions == "><") && dist <= args.fold {
        if args.no_sym {
              break_counts[0] +=1; 
            return String::from("Foldback"); 
        } else {
            return check_sym(brk, args, break_counts); 
        }
    } 

    break_counts[2] +=1; 
    return String::from("Pass"); 
}

// fucntion to check whether the foldback artifact occurs nearly in the middle of the read 
fn check_sym(brk:&Breakpoint, args: &Cli, break_counts: &mut Vec<u64>) -> String {
    //consider symetric read if break occurs +/- 5% of middle of read 
    let rlen = brk.read_len as f32; 
    let r1 = brk.r2_loc as f32;
    let r2 = brk.r2_loc as f32; 
     
    //consider middle of r1 and r2 coords to be the break location in read coordinates
    let r_break_ave = (r1 + r2) / 2.0; 

    let range_min = rlen/2.0 - (args.margin * rlen); 
    let range_max = rlen/2.0 + (args.margin * rlen); 
    if range_min <= r_break_ave  && r_break_ave <= range_max {
        break_counts[0] +=1; 
        return String::from("Foldback");
    } else {
        break_counts[2] +=1; 
        return String::from("Pass"); 
    }
} 

//get the read level classification based on the classification of the breakpoints in the read and update all counts 
pub fn read_level_class(label_counts: &Vec<u32>, cur_counts: &mut Vec<u64>) {
    //if only pass breaks found, read is not artifact
    if label_counts[2] > 0 && (label_counts[0] + label_counts[1] == 0){
        cur_counts[2] +=1; 
    } else if label_counts[0] >=  label_counts[1]{ 
        //more folds breaks
        //if even one fold or chim found, read is artifact and classified by 
        //which one was more frequent, a tie goes to fold
        cur_counts[0] +=1; 
    }else {
        //more chim breaks
        cur_counts[1] +=1; 
    }
}


//print summary statistics to stdout
pub fn print_report(reads_pass_filter: u64, read_counts: &Vec<u64>, break_counts: &Vec<u64>, args: &Cli ) -> io::Result<()>   {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    writeln!(handle, "{}", "*".repeat(100))?;
    writeln!(handle, "Breakinator summary report:")?;
    writeln!(handle, "Command: {}", env::args().collect::<Vec<_>>().join(" "))?;
    writeln!(handle, "Filtering Criteria: MapQ \u{2265} {} and min_alignment_len \u{2265} {}", args.min_mapq, args.min_map_len)?;
    writeln!(handle, "\nResults\n{}", "-".repeat(20))? ;
    writeln!(handle, "Num reads passed filter: {}", add_commas(reads_pass_filter))? ;
    writeln!(handle, "Num breakpoints detected: {} on {} unique reads" , add_commas(break_counts.iter().sum::<u64>()), add_commas(read_counts.iter().sum::<u64>()))? ;
    writeln!(handle, "\nFoldback artifacts:")? ; 
    writeln!(handle, "Num Foldback READS detected: {}  ({}% of all reads)" , add_commas(read_counts[0]), get_percent(read_counts[0], reads_pass_filter))? ;
    writeln!(handle, "Num Foldback BREAKPOINTS detected: {}  ({}% of all breakpoints)" , add_commas(break_counts[0]), get_percent(break_counts[0] ,  break_counts.iter().sum::<u64>()))? ;
    writeln!(handle, "\nChimeric artifacts:")? ; 
    writeln!(handle, "Num Chimeric READS detected: {}  ({}% of all reads)" , add_commas(read_counts[1]), get_percent(read_counts[1], reads_pass_filter))? ;
    writeln!(handle, "Num chimeric BREAKPOINTS detected: {}  ({}% of all breakpoints)" , add_commas(break_counts[1]), get_percent(break_counts[1], break_counts.iter().sum::<u64>()))? ;
    writeln!(handle, "{}", "*".repeat(100))?;
    return Ok(())
}

//print summary statistics to stdout in table format
pub fn print_table(reads_pass_filter: u64, read_counts: &Vec<u64>, break_counts: &Vec<u64>, file_name: String ) -> io::Result<()>   {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let tot_breaks = break_counts.iter().sum::<u64>(); 
    let all_stats = vec![reads_pass_filter.to_string(), tot_breaks.to_string(), read_counts.iter().sum::<u64>().to_string(), read_counts[0].to_string(), get_percent(read_counts[0],reads_pass_filter) + "%", break_counts[0].to_string(), get_percent(break_counts[0],tot_breaks)+ "%", read_counts[1].to_string(), get_percent(read_counts[1],reads_pass_filter)+ "%", break_counts[1].to_string(), get_percent(break_counts[1],tot_breaks) + "%", file_name]; 
    writeln!(handle,"#Reads_passed\tall_break\tUniq_artifact_reads\tFold_reads\tFold_reads%\tFold_breaks\tFold_breaks%\tChim_reads\tChim_reads%\tChim_breaks\tChim_breaks%\tsample")?; 
    writeln!(handle, "{}", all_stats.join("\t"))?; 
    Ok(())
}

// convert a number to a comma seperated String
fn add_commas(num: u64) -> String {
    return num.to_string()
    .as_bytes()
    .rchunks(3)
    .rev()
    .map(std::str::from_utf8)
    .collect::<Result<Vec<&str>, _>>()
    .unwrap()
    .join(","); 
}

// calculate the precentage and return as String 
fn get_percent(a:u64, b:u64) -> String {
    if b == 0 {
        return String::from("n/a")
    }
    let p = 100.0 * a as f32 / b as f32; 
    return format!("{:.3}", p); 
}
