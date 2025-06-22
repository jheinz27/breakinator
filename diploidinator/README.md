# The diploid-inator

Most aligners were not designed for diploid assemblies(eg. [HG002](https://github.com/marbl/HG002)), so when aligning reads to a diploid assembly, the mapping quality for reads may be lower, as there are multiple locations the read can align to well. We have developed a simple script to align reads to each haploid of the diploid assembly and thne parse both paf files to choose the better alignment of the read based on the alignment score

---

## Installation

(Recommended) You can download the latest compiled binary from the [Releases](https://github.com/jheinz27/diploidinator/releases) page.
```
wget https://github.com/jheinz27/diploidinator/releases/download/v0.1.0/diploidinator
chmod +x diploidinator
./diploidinator
```

or build from source 
``` 
git clone https://github.com/jheinz27/diploidinator.git
cd diploidinator
cargo build --release
./target/release/diploidinator
```

## Example Usage
NOTE: It is important to use the `--secondary=no --paf-no-hit` flags when aligning with [minimap2](https://github.com/lh3/minimap2). The diploidinator currently only works on paf files. 
```
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.MATERNAL.fasta read.fastq > out_mat.paf
minimap2 -cx splice -uf -k14 -t 16 --secondary=no --paf-no-hit hg002v1.1.PATERNAL.fasta reads.fastq > out_pat.paf 
diploidinator out_mat.paf out_pat.paf > out_haps_merge.paf
```

