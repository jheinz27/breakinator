# The Breakinator

<img src="https://github.com/user-attachments/assets/447a923e-c4d1-4331-8a81-130f48144ca0" alt="Breakinator Logo" width="300" align="left" hspace="20"/>

### Table of Contents
[Installation](#installation)<br>
&emsp;[Prebuilt Binaries](#prebuilt-binaries)<br>
&emsp;[Bioconda](#bioconda)<br>
&emsp;[Install from source](#install-from-source)<br>
[Breakinator Usage](#breakinator-usage)<br>
&emsp;[Example Usage](#example-usage)<br>
&emsp;&emsp;[For SAM/BAM/CRAM](#for-sambamcram)<br>
&emsp;&emsp;[Coordinate sorted input files (`-s`)](#coordinate-sorted-input-files--s)<br>
&emsp;&emsp;[For PAF](#for-paf-include---paf-flag)<br>
[Generating PAF files](#generating-paf-files)<br>
[Optional: turn off symmetry filter for foldback artifacts](#optional-turn-off-symmetry-filter-for-foldback-artifacts)<br>
[Merging Breakpoints Into Consensus Locations](#merging-breakpoints-into-consensus-locations)<br>
[Citation](#citation)

<br clear="left"/>

The Breakinator identifies and flags putative artifact reads (foldbacks and chimeric) by parsing SAM/BAM/CRAM or PAF alignment files.

## Installation

### Prebuilt Binaries

```
wget https://github.com/jheinz27/breakinator/releases/download/v1.2.0/breakinator-v1.2.0-x86_64-unknown-linux-musl.tar.gz
tar -xvzf breakinator-v1.2.0-x86_64-unknown-linux-musl.tar.gz
./breakinator-v1.2.0-x86_64-unknown-linux-musl/bin/breakinator --help
```

Other versions, including macOS, are available on the [Releases](../../releases) page.

### Bioconda
```
conda install -c bioconda -c conda-forge breakinator
```

### Install from source
```
git clone https://github.com/jheinz27/breakinator
cd breakinator/breakinator
cargo build --release
./target/release/breakinator --help
```


#### Prerequisites
- Rust programming language >= v1.70
- clap = "4.0" 
- rust-htslib = "0.46.0"
## Breakinator Usage
```
Flag foldbacks and chimeric reads from SAM/BAM/CRAM or PAF input

Usage: breakinator [OPTIONS] --input <FILE>

Options:
  -i, --input <FILE>       Input SAM/BAM/CRAM/PAF file (grouped by read ID unless -s is given)
  -s, --sorted             Input is coordinate sorted
      --paf                Input file is PAF
  -q, --min-mapq <INT>     Minimum mapping quality [default: 10]
  -a, --min-map-len <INT>  Minimum alignment length (bps) [default: 200]
      --no-sym             Report all foldback reads, not just those with breakpoint within margin of middle of read
  -g, --genome <FASTA>     Reference genome FASTA used (must be provided for CRAM input)
  -m, --margin <FLOAT>     [0-1], Proportion from center of read on either side to be considered sym foldback artifact [default: 0.1]
      --rcoord             Print read coordinates of breakpoint in output
  -o, --out <FILE>         Output file name [default: breakinator_out.txt]
  -c, --chim <INT>         Minimum distance to be considered chimeric [default: 1000000]
  -f, --fold <INT>         Max distance to be considered foldback [default: 200]
      --tabular            Print a TSV table instead of the default report (useful if evaluating multiple samples)
  -t, --threads <INT>      Number of threads to use for BAM/CRAM I/O [default: 2]
  -h, --help               Print help
  -V, --version            Print version
```
### Example Usage 

By default, Breakinator expects alignments grouped by read ID (the default output of minimap2), and processes one read at a time to keep memory low. For coordinate-sorted SAM/BAM/CRAM files, use `-s` (see below).

#### For SAM/BAM/CRAM
```
minimap2 -ax map-ont  genome.fa reads.fastq > alignments.sam
./breakinator -i alignments.sam -o breakinator_out.txt
```
Input can also be streamed from stdin with `-i -`, e.g. directly from minimap2 without writing a SAM file:
```
minimap2 -ax map-ont genome.fa reads.fastq | ./breakinator -i - -o breakinator_out.txt
```

#### Coordinate sorted input files (`-s`)
As of version 1.2 Breakinator supports coordinate sorted input files with the `--sorted / -s` flag (requires additional memory, typically less than 100 MB extra)
```
./breakinator -i alignments_sorted.bam -s -o breakinator_out.txt
./breakinator -i alignments_sorted.cram -s -g genome.fa -o breakinator_out.txt
```

#### For PAF (include `--paf` flag)
```
minimap2 -cx map-ont --secondary=no genome.fa reads.fastq > alignments.paf
./breakinator -i alignments.paf --paf -o breakinator_out.txt
```
PAF input must still be grouped by read name. If PAF has been sorted by position, we recommend sorting it by read name and piping the result directly:
```
sort -k1,1 input.paf | ./breakinator --paf -i - -o breakinator_out.txt
```

## Generating PAF files

The Breakinator can also handle PAF files to input to the Breakinator. To generate these, we recommend using [minimap2](https://github.com/lh3/minimap2) with the `-c` and `--secondary=no` parameters. Secondary alignments will be ignored by the Breakinator, however including them will increase the processing time. 

Example:
```
minimap2 -cx map-ont --secondary=no genome.fa reads.fastq > alignments.paf
```

The PAF can also be generated by converting a SAM file to a PAF with [paftools.js](https://github.com/lh3/minimap2/blob/master/misc/paftools.js) using the `-p` parameter. 

Example: 
```
paftools.js sam2paf -p alignments.sam > alignments.paf
```


## Optional: turn off symmetry filter for foldback artifacts

If running on a sample where you want to investiage all potential foldback events, we recommend turning off the symmetry filter with the `--no-sym` flag. 

```
./breakinator -i alignments.paf --no-sym
```
<img width="742" alt="Screenshot 2025-05-09 at 10 15 35 AM" src="https://github.com/user-attachments/assets/c66855bb-5fbd-4143-a884-9bd200a4395f" />

## Preprocessing for alignment to diploid  assemblies with Diploidinator(DEPRECATED)

**Diploidinator has been deprecated** 

We recommend using [HipHap](https://github.com/jheinz27/hiphap) for diploid genome alignments.


## Merging Breakpoints Into Consensus Locations
To evaluate how many unique breakpoints are in the sample and how much read support they have, we developed a simple script to merge breakpoints together if they occur within 100bps (default `-w`) of eachother. We require at least 2 reads (default `-s`) of support to report a consensus breakpoint location. 
```
usage: merge_breaks.py [-h] -i <breakpoints.txt> [-w --merge_window] [-s --min_support] > merged_breaks.txt

Merge Break Points from Breakinator output

optional arguments:
  -h, --help            show this help message and exit
  -i <breakpoints.txt>  input breakinator stdout
  -w --merge_window     Size of window to merge break points in
  -s --min_support      minimum reads supporting breakpoint
```

## Citation
If the Breakinator has helped you in your research, please cite our [paper](https://link.springer.com/article/10.1186/s12864-025-12492-y).

Heinz, J.M., Meyerson, M. & Li, H. "Detecting foldback artifacts in long-reads" *BMC Genomics* 27, 144 (2026). doi:[10.1186/s12864-025-12492-y](https://doi.org/10.1186/s12864-025-12492-y)


