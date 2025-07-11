import sys
import numpy as np 
import argparse

#funciton that takes a cluster of breakpoints within the margin to merge into one line
def get_merges(cluster,min_support): 
    if len(cluster) < min_support: 
        return None 
    #convert clsuer into array
    arr = np.array(cluster)
    #take break point locations to be the median of the locations in the cluster of reads
    start_val = round(np.median([int(i) for i in arr[:,1]]))
    end_val = round(np.median([int(i) for i in arr[:,4]])) 
    return '\t'.join([arr[0,0], str(start_val), arr[0,2], arr[0,3], str(end_val), ','.join(arr[:,5]), ','.join(arr[:,6]), str(len(arr[:,6]))]) 
   
#convert the breakpoint text format into vcf format 
def vcf_format(original, number): 
    out = [] 
    s = original.strip().split('\t')
    qual = str(round(np.mean([int(x) for x in s[5].split(',')])))
    read_count = str(len(s[5].split(';')))
    info = ';READ_COUNTS='+ read_count + ';READ_IDS=' + s[6]
    if s[2] == '>>':
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', 'N[' + s[3]+':'+ s[4]+'[', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1)  + info]))
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', ']' + s[0]+':'+ s[1]+']N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number)  + info]))
    elif s[2] == '><':
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', 'N]' + s[3]+':'+ s[4]+']', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1) + info])) 
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', 'N]' + s[0]+':'+ s[1]+']', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number) + info]))
    else: 
        out.append('\t'.join([s[0],s[1],'bnd_' + str(number), 'N', '[' + s[3]+':'+ s[4]+'[N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number+1) + info]))
        out.append('\t'.join([s[3],s[4],'bnd_' + str(number+1), 'N', '[' + s[0]+':'+ s[1]+'[N', qual, '.', 'SVTYPE=BND;MATEID=bnd_'+ str(number) + info]))

    return out

def identify_diff_breaks(group, margin): 
    sorted_group = sorted(group, key = lambda x: (int(x[4])))
    i = 0 
    groups = [] 
    while i < len(sorted_group): 
        j = 0 
        while i+j+1 < len(sorted_group) and abs(int(sorted_group[i+j][4]) - int(sorted_group[i+j+1][4])) <= margin: 
            j +=1
        groups.append(sorted_group[i:i+j+1]) 
        i+= (j+1)
    return groups 

def merge_breaks(breakpoints, margin=100, support=2): 
    all_merges= []   
    sorted_lines = sorted(breakpoints, key = lambda x: (x[0],x[3], int(x[1]))) 
    i = 0 
    #group breakpoints to merge 
    while i < len(sorted_lines): 
        j = 0 
        while i+j+1 < len(sorted_lines) and sorted_lines[i+j][0] == sorted_lines[i+j+1][0] and abs(int(sorted_lines[i+j][1]) - int(sorted_lines[i+j+1][1])) <= margin and sorted_lines[i+j][3] == sorted_lines[i+j+1][3]: 
            j +=1
            
        groupings = identify_diff_breaks(sorted_lines[i:i+j+1], margin) 
        for group in groupings:
            out = get_merges(group, support)
            if out:
                all_merges.append(out)
        i += (j + 1)  
    return all_merges

if __name__ == "__main__": 
    parser = argparse.ArgumentParser(description='Merge Break Points from Breakinator output')
    parser.add_argument('-i', metavar='<breakpoints.txt>', required=True, help='input breakinator stdout')
    parser.add_argument('-w', metavar= '--merge_window', required= False, type=int, default = 100, help = 'Size of window to merge break points in')
    parser.add_argument('-s', metavar= '--min_support', required= False, type=int, default = 2, help = 'minimum reads supporting breakpoint') 
   
    args = parser.parse_args()
    brkFile = args.i
    margin = args.w
    min_support = args.s
   
    lines = [] 
    with open(brkFile, 'r') as f:
        for line in f:
            lines.append(line.strip().split('\t'))
   
    sys.stdout.write('\n'.join(merge_breaks(lines, margin, min_support))) 



