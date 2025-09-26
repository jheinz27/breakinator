use std::{env, iter::Peekable, cmp::min, io::{self, Write}};
pub mod cli;
pub use cli::Cli;


//function to move ahead one read group at a time and update cluster with the group
pub fn get_clusters<I>(lines: &mut Peekable<I>, cluster: &mut Vec<String>)-> io::Result<()>
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
    return Ok(())
}

//filter out read alignments that fail length or mapQ filters or is secondary alignment 
pub fn filter_alignments(all_maps: &Vec<String>, mapq: i32, map_len: i32) -> Vec<Vec<&str>>{
    let mut passed_filter: Vec<Vec<&str>> =  Vec::new(); 
    for map in all_maps{
        let fields: Vec<&str> = map.split('\t').collect();
        if fields[11].parse::<i32>().expect("MAPQ field was not a valid integer") >= mapq && (fields[3].parse::<i32>().expect("Alignment end not a valid integer")- fields[2].parse::<i32>().expect("Alignment start not a valid integer")) >= map_len {
            passed_filter.push(fields);
        }
    }
    return passed_filter; 
}

//find the breakpoints and a read that has at least one supplementary alignment
pub fn determine_break(clust: &mut Vec<Vec<&str>>, args: &Cli, read_counts: &mut Vec<i64>, break_counts: &mut Vec<i64>) ->  Result<Vec<Vec<String>>,  String> {
    //sort by start location of aligment in read 
    clust.sort_by_key(|line| {line[2].parse::<i32>().expect("invalid int in field 3")});
    let mut out: Vec<Vec<String>> = Vec::new(); 
    let mut labels: Vec<i32> = vec![0; 3];  //[fold, chim, pass count]

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
            _ => {return Err(format!( "determine_break: unexpected strand: {}",clust[i].join("\t")))}
        };
    
        let b2 = match clust[i+1][4] {
            //second break on forward strand, take start of later alignment 
            "+" => {directions.push('>'); vec![clust[i+1][5], clust[i+1][7]]} 
            //second break on reverse strand, take end of later alignment
            "-" => {directions.push('<'); vec![clust[i+1][5], clust[i+1][8]]}
            _ => {return Err(format!( "determine_break: unexpected strand: {}",clust[i].join("\t")))}
        }; 
        
        //set mapq val to be min of mapq value for both sides of breakpoint
        let mapq = min(clust[i][11].parse::<i32>().expect("invalid int in field 11"), clust[i+1][11].parse::<i32>().expect("invalid int in field 11"));
        let mut break_info = vec![b1[0].to_string(), b1[1].to_string(), directions, b2[0].to_string(), b2[1].to_string(),mapq.to_string(), clust[0][0].to_string()];  

        let rlen = clust[i][1].parse::<f32>().expect("invalid int in field 2"); 
        let rbreak = clust[i][3].parse::<f32>().expect("invalid int in field 4"); 

        //get artifact or pass classification 
        let label = classify_break( &break_info, args, rlen, rbreak, break_counts); 
       
        match label.as_str() {
            "Foldback" => labels[0] += 1,
            "Chimeric" => labels[1] += 1,
            _ => labels[2] += 1,
        }
        break_info.push(label); 

        //Add read coordinates of breakpoint
        if args.rcoord {
            break_info.push(clust[i][3].to_string()); 
            break_info.push(clust[i+1][2].to_string()); 
        }
    
        out.push(break_info); 
    
    }
    //classify the read based on the classification of all breakpoints on the read
    read_level_class(&labels, read_counts); 
    return Ok(out); 
    
} 

//print summary statistics to stdout
pub fn print_report(reads_pass_filter: i64, read_counts: &Vec<i64>, break_counts: &Vec<i64>, args: &Cli ) -> io::Result<()>   {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    writeln!(handle, "{}", "*".repeat(100))?;
    writeln!(handle, "Breakinator summary report:")?;
    writeln!(handle, "Command: {}", env::args().collect::<Vec<_>>().join(" "))?;
    writeln!(handle, "Filtering Criteria: MapQ \u{2265} {} and min_alignment_len \u{2265} {}", args.min_mapq, args.min_map_len)?;
    writeln!(handle, "\nResults\n{}", "-".repeat(20))? ;
    writeln!(handle, "Num reads passed filter: {}", add_commas(reads_pass_filter))? ;
    writeln!(handle, "Num breakpoints detected: {} on {} unique reads" , add_commas(break_counts.iter().sum::<i64>()), add_commas(read_counts.iter().sum::<i64>()))? ;
    writeln!(handle, "\nFoldback artifacts:")? ; 
    writeln!(handle, "Num Foldback READS detected: {}  ({}% of all reads)" , add_commas(read_counts[0]), get_percent(read_counts[0], reads_pass_filter))? ;
    writeln!(handle, "Num Foldback BREAKPOINTS detected: {}  ({}% of all breakpoints)" , add_commas(break_counts[0]), get_percent(break_counts[0] ,  break_counts.iter().sum::<i64>()))? ;
    writeln!(handle, "\nChimeric artifacts:")? ; 
    writeln!(handle, "Num Chimeric READS detected: {}  ({}% of all reads)" , add_commas(read_counts[1]), get_percent(read_counts[1], reads_pass_filter))? ;
    writeln!(handle, "Num chimeric BREAKPOINTS detected: {}  ({}% of all breakpoints)" , add_commas(break_counts[1]), get_percent(break_counts[1], break_counts.iter().sum::<i64>()))? ;
    writeln!(handle, "{}", "*".repeat(100))?;
    return Ok(())
}

//print summary statistics to stdout in table format
pub fn print_table(reads_pass_filter: i64, read_counts: &Vec<i64>, break_counts: &Vec<i64>, file_name: String ) -> io::Result<()>   {
    let stdout = io::stdout();
    let mut handle = stdout.lock();
    let tot_breaks = break_counts.iter().sum::<i64>(); 
    let all_stats = vec![reads_pass_filter.to_string(), tot_breaks.to_string(), read_counts.iter().sum::<i64>().to_string(), read_counts[0].to_string(), get_percent(read_counts[0],reads_pass_filter) + "%", break_counts[0].to_string(), get_percent(break_counts[0],tot_breaks)+ "%", read_counts[1].to_string(), get_percent(read_counts[1],reads_pass_filter)+ "%", break_counts[1].to_string(), get_percent(break_counts[1],tot_breaks) + "%", file_name]; 
    writeln!(handle,"#Reads_passed\tall_break\tUniq_artifact_reads\tFold_reads\tFold_reads%\tFold_breaks\tFold_breaks%\tChim_reads\tChim_reads%\tChim_breaks\tChim_breaks%\tsample")?; 
    writeln!(handle, "{}", all_stats.join("\t"))?; 
    Ok(())
}

// convert a number to a comma seperated String
fn add_commas(num: i64) -> String {
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
fn get_percent(a:i64, b:i64) -> String {
    if b == 0 {
        return String::from("n/a")
    }
    let p = 100.0 * a as f32 / b as f32; 
    return format!("{:.3}", p); 
}

// fucntion to check whether the foldback artifact occurs nearly in the middle of the read 
fn check_sym(read_len: f32, read_break:f32, args: &Cli, break_counts: &mut Vec<i64>) -> String {
    //consider symetric read if break occurs +/- 5% of middle of read 
    let range_min = read_len/2.0 - (args.margin * read_len); 
    let range_max = read_len/2.0 + (args.margin * read_len); 
    if range_min <= read_break  && read_break <= range_max {
        break_counts[0] +=1; 
        return String::from("Foldback");
    } else {
        break_counts[2] +=1; 
        return String::from("Pass"); 
    }
} 

//classify break as either chimeric, foldback, or pass
fn classify_break(brk: &Vec<String>,  args: &Cli, read_len: f32, read_break: f32, break_counts: &mut Vec<i64>) -> String {
    let dist = (brk[4].parse::<i32>().expect("cannot be parsed to int") - brk[1].parse::<i32>().expect("cannot be parsed to int")).abs(); 
    if brk[0] != brk[3] {
        break_counts[1] +=1; 
        return String::from("Chimeric"); 
    } else if dist >= args.chim {
        break_counts[1] +=1; 
        return String::from("Chimeric"); 
    } else if (brk[2] == "<>" || brk[2] == "><") && dist <= args.fold {
        if args.sym {
            return check_sym(read_len, read_break, args, break_counts); 
        } else {
            break_counts[0] +=1; 
            return String::from("Foldback"); 
        }
    } 
    break_counts[2] +=1; 
    return String::from("Pass"); 
}

//get the read level classification based on the classification of the breakpoints in the read and update all counts 
fn read_level_class(label_counts: &Vec<i32>, cur_counts: &mut Vec<i64>) {
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
