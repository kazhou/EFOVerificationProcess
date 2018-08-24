#!/bin/sh
# Grid Engine options
#$ -N  matchRel
#$ -cwd
#$ -S /bin/bash
#$ -V

module load sge
module load R/3.5.0

read -p 'Enter target(s) to process ("all" for all targets; "approved" for approved targets only;"phase3" for phase 3 & beyond targets;
    <filename.txt> for targets from file; indiv. gene symbols sep. with comma, no space, upper case;): ' targets
read -p 'Enter threshold: ' threshold
read -p 'Enter desired P-value: ' pval
read -p 'Normal p-value or FDR-adjusted p-value? ("normal" or "FDR"): ' pstat
read -p 'Enter number of rows to show ("all" for all): ' nrows
read -p 'Enter date FET files were generated ("YYYYMMDD"): ' date

Rscript ~/scripts/03_matchWrapper.R $targets $threshold $pval $pstat $nrows $date