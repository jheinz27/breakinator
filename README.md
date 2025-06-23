
<div style="text-align: center;">
  <img src="https://github.com/user-attachments/assets/447a923e-c4d1-4331-8a81-130f48144ca0" alt="Breakinator Logo" width="300"/>
</div>

# The Breakinator
The Breakinator identifies and flags putative artifact reads (foldbacks and chimeric) by parsing a PAF file.

## Installation

```
git clone git@github.com:jheinz27/breakinator.git
cd breakinator
python breakinator.py -h
```

#### Prerequisites

- Python 3.7 or higher  


## Usage
```
usage: breakinator.py [-h] -i FILE [-m INT] [-a INT] [--sym] [--margin FLOAT] [-o FILE] [--chim INT] [--fold INT] [--tabular]

Flag foldbacks and chimeric reads from PAF input

optional arguments:
  -h, --help      show this help message and exit
  -i FILE         PAF file sorted by read IDs
  -m INT          Minimum mapping quality (integer)
  -a INT          Minimum alignment length (bps)
  --sym           Only report palindromic foldback reads within margin
  --rcoord        Print read coordinates of breakpoint in output
  --margin FLOAT  [0-1], With --sym, Proportion from center on either side to be considered foldback artifact
  -o FILE         Output file name
  --chim INT      Minimum distance to be considered chimeric
  --fold INT      Max distance to be considered foldback
  --tabular       Return report as a tsv file (useful for evaluating multiple files)
```


## Optional symmetry filter for foldback artifacts

If running on a sample where you expect there may be true foldback events (e.g cancer data) we recommend using the `--sym` flag to only consider reads that have the foldback withing 5% (change with `--margin`) on either side of the middle of the read as foldback artifacts. See diagram below. 

```
python breakinator.py -i alignments.paf --sym --margin 0.03
```
<img width="742" alt="Screenshot 2025-05-09 at 10 15 35 AM" src="https://github.com/user-attachments/assets/c66855bb-5fbd-4143-a884-9bd200a4395f" />

## Citation
If this tool has helped you in your research, please cite our preprint at: 


