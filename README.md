<div style="text-align: center;">
  <img src="https://github.com/user-attachments/assets/447a923e-c4d1-4331-8a81-130f48144ca0" alt="Breakinator Logo" width="400"/>
</div>

The Breakinator identifies and flags putative artifact reads (foldbacks and chimeric) by parsing a PAF file.

##installation


## usage
usage: breakinator.py [-h] -i FILE [-m INT] [-a INT] [--sym] [--margin FLOAT] [-o FILE] [--chim INT] [--fold INT] [--tabular]

Flag foldbacks and chimeric reads from PAF input

optional arguments:
  -h, --help      show this help message and exit
  -i FILE         PAF file sorted by read IDs
  -m INT          Minimum mapping quality (integer)
  -a INT          Minimum alignment length (bps)
  --sym           Only report palindromic foldback reads within margin
  --margin FLOAT  [0-1], With --sym, Proportion from center on either side to be considered foldback artifact
  -o FILE         Output file name
  --chim INT      Minimum distance to be considered chimeric
  --fold INT      Max distance to be considered foldback
  --tabular       Return report as a tsv file (useful for evaluating multiple files)
