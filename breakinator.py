import sys 
import argparse

#function to check if foldback occurs approximately symetrically in the read, if not, label read as Pass
def check_sym(read_len, read_break): 
    #consider symetric read if break occurs +/- 5% of middle of read
    read_len_range = [read_len/2 - fold_margin * read_len, read_len/2 + fold_margin * read_len]
    if read_len_range[0] < read_break < read_len_range[1]: 
        break_counts[1] += 1 
        return 'Foldback'
    else: 
        return 'Pass'

#function to detect and return label of technical aritfact
def check_artifact(brk,read_len, read_break):
    break_counts[0] += 1 
    if brk[0] != brk[3]: 
        break_counts[2] += 1
        return "Chimeric" 
    elif abs(int(brk[4]) - int(brk[1])) >= min_chim: 
        break_counts[2] += 1 
        return "Chimeric" 
    elif (brk[2] == '<>' or brk[2] == '><') and abs(int(brk[4]) - int(brk[1])) <= max_fold:
        #check if entire read is approximately palindromic if sym filter specified 
        if sym_filter: 
            return check_sym(read_len, read_break)
        else:
            break_counts[1] += 1 
            return 'Foldback'
    else: 
        return 'Pass'

#function to take in all split alignments of one read
#returns break point in read 
def breakpoint(cluster): 
    #sort by start location of aligment in read 
    sort = sorted(cluster, key = lambda x: int(x[2]))
    all_labels = []
    outs = []
    #check every concurrent alignment in a read with n split alignments
    for i in range(len(sort)-1): 
        #get strands of both aligments
        strand = [sort[i][4], sort[i+1][4]]
        directions = [] 
        #First break on forward strand, take end of first alignment
        if strand[0] == '+':
            directions.append('>')
            b1 = [sort[i][5],sort[i][8]] 
        #First break on reverse strand, take start of first alignment
        else:
            directions.append('<')
            b1 = [sort[i][5],sort[i][7]]
        #second break on forward strand, take start of later alignment 
        if strand[1] == '+': 
            directions.append('>')
            b2 = [sort[i+1][5],sort[i+1][7]] 
        #second break on reverse strand, take end of later alignment 
        else: 
            directions.append('<')
            b2 = [sort[i+1][5],sort[i+1][8]] 
        directs = ''.join(directions)

        brk_info = [b1[0], b1[1], directs, b2[0], b2[1],str(min([int(sort[i][11]),int(sort[i+1][11])])),cluster[0][0]]
        label = check_artifact(brk_info, int(sort[i][1]), int(sort[i][3])) 
        all_labels.append(label)
        brk_info.append(label)
        if rcoord: 
            brk_info.append(sort[i][3])
            brk_info.append(sort[i+1][2])
        outs.append(brk_info)

    return outs, all_labels

def get_percent(a,b): 
    return str(round(100*a/b,3)) + '%'

#function to generate terminal output from 
def make_report(stats):
    report = ['*'*100,'Summary report:']
    report.append('Command: ' + ' '.join(sys.argv))
    report.append(f'Filtering criteria: MapQ >= {min_mapQ} and min_alignment_len >{min_map_len}\n\nResults:\n' + '-' *15)
    report.append(f'Num reads passed filter: {stats[0]:,}') 
    report.append(f'Num breakpoints detected: {stats[1]:,} on {stats[2]:,} unique reads')
    report.append('\nFoldback artifacts:')
    report.append(f'Num Foldback READS detected: {stats[3]:,} ({stats[4]} of all reads)')
    report.append(f'Num Foldback BREAKPOINTS detected: {stats[5]:,} ({stats[6]} of all breakpoints)' )
    report.append('\nChimeric artifacts:')
    report.append(f'Num Chimeric READS detected: {stats[7]:,} ({stats[8]} of all reads)')
    report.append(f'Num Chimeric BREAKPOINTS detected: {stats[9]:,} ({stats[10]} of all breakpoints)') 
    report.append('*'*100)
    return report 
 
#function to assign read artifact read label if multiple breaks on read
def update_read_labels(labels): 
    labels_counts= [0,0,0]
    #get counts for every breakpoint in read
    for lab in labels: 
        if lab == 'Foldback': 
            labels_counts[0] +=1 
        elif lab == 'Chimeric': 
            labels_counts[1] +=1 
        else: 
            labels_counts[2] +=1
    #if only pass breaks found, read is not artefact
    if labels_counts[2] > 0 and (labels_counts[0] + labels_counts[1] == 0):
        return  
    #if even one fold or chim found, read is artifact and classified by 
    #which one was more frequent, a tie goes to fold
    #more folds breaks
    if labels_counts[0] >=  labels_counts[1]: 
        read_counts[1] +=1  
    #more chim breaks
    elif labels_counts[0] < labels_counts[1]: 
        read_counts[2] +=1  
    return 

def main(paf):
    passed_filter_ids = []
    #read input PAF
    mapped_read_count = 0 
    all_breaks = []
    cur = ''
    clust = [] 
    with open(paf, 'r') as f: 
        lowQ = True
        for line in f:
            lowQ = False
            s = line.strip().split('\t')
            #skip read alignments that fail length or mapQ filters or is secondary alignment 
            if int(s[11]) < min_mapQ or (int(s[3])  - int(s[2]) <= min_map_len) :
                lowQ = True
                continue
            #group aligments of read by read ID
            if s[0].strip() == cur: 
                clust.append(s)
            #once all alignments are gathered evaluate the breakpoint
            else:
                mapped_read_count += 1 
                passed_filter_ids.append(s[0].strip())
                if len(clust) > 1: 
                    out, labels = breakpoint(clust)
                    read_counts[0] += 1 
                    update_read_labels(labels) 
                    for o in out:
                        all_breaks.append(o)
                clust = [s]
            cur = s[0] 
    #check last cluster of read aligments if exists
    if not lowQ: 
        mapped_read_count += 1
        passed_filter_ids.append(s[0].strip())
    if len(clust) > 1:
        out, labels = breakpoint(clust)
        read_counts[0] += 1 
        update_read_labels(labels)
        for o in out:
            all_breaks.append(o)

    #write read ID and break point, and artifact flag
    with open(outFile, 'w') as read_results:
        for b in all_breaks: 
            read_results.write('\t'.join(b) + '\n')

    #add in a check for when reads passing our filter is 0 ! 
    
    #[#Reads_passed,all_break,Uniq_artifact_reads,Fold_reads, Fold_reads%, Fold_breaks,Fold_breaks%,Chim_reads,Chim_reads%,Chim_breaks,Chim_breaks%, sample]'
    vals = [mapped_read_count, break_counts[0],read_counts[0], read_counts[1], get_percent(read_counts[1], mapped_read_count),  break_counts[1],get_percent(break_counts[1], break_counts[0]),
             read_counts[2], get_percent(read_counts[2], mapped_read_count), break_counts[2], get_percent(break_counts[2], break_counts[0]), pafFile]
    if tabular:
        sys.stdout.write('#Reads_passed\tall_break\tUniq_artifact_reads\tFold_reads\tFold_reads%\tFold_breaks\tFold_breaks%\tChim_reads\tChim_reads%\tChim_breaks\tChim_breaks%\tsample\n')
        sys.stdout.write('\t'.join(str(v) for v  in vals) + '\n')  
    else:
        report_out = make_report(vals)
        sys.stdout.write('\n'.join(report_out) + '\n')

    with open('reads_passed_our_filter.txt', 'w') as o: 
        o.write('\n'.join(passed_filter_ids))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Flag foldbacks and chimeric reads from PAF input')
    parser.add_argument('-i', metavar='FILE', required=True, help='PAF file sorted by read IDs')
    parser.add_argument('-m', metavar='INT', required= False, type=int, default = 10, help = 'Minimum mapping quality (integer)')
    parser.add_argument('-a', metavar='INT', required= False, type=int, default = 200, help = 'Minimum alignment length (bps)')
    parser.add_argument('--sym', action= "store_true",  help= 'Only report palindromic foldback reads within margin') 
    parser.add_argument('--rcoord', action= "store_true",  help= 'Print read coordinates of breakpoint in output') 
    parser.add_argument('--margin', metavar='FLOAT', required = False, type=float, default = 0.05, help='[0-1], With --sym, Proportion from center on either side to be considered foldback artifact' ) 
    parser.add_argument('-o', metavar='FILE', required= False, type=str, default ='stdout.txt', help= 'Output file name' )
    parser.add_argument('--chim', metavar='INT', required= False, type=int, default = 1_000_000, help = 'Minimum distance to be considered chimeric')
    parser.add_argument('--fold', metavar='INT', required= False, type=int, default = 200, help = 'Max distance to be considered foldback')
    parser.add_argument('--tabular', action= "store_true", help= 'Return report as a tsv file (useful for evaluating multiple files)')

    args = parser.parse_args()
    pafFile = args.i
    min_mapQ = args.m
    min_map_len = args.a
    sym_filter = args.sym 
    outFile = args.o
    min_chim = args.chim
    max_fold = args.fold
    fold_margin = args.margin
    tabular = args.tabular
    rcoord = args.rcoord

    if not 0<= fold_margin <=1: 
        parser.error('-margin must be [0-1]')
    
    global read_counts
    global break_counts
    read_counts = [0,0,0] #[tot_reads_break, tot_reads_fold, tot_reads_chim]
    break_counts = [0,0,0] #[tot_breaks, tot_folds, tot_chimeras]
    main(pafFile)
