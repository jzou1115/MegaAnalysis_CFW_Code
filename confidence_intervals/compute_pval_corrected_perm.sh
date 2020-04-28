#!/bin/sh
#$ -N perm
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=10G,highp,h_rt=23:59:59
#$ -m a
#$ -M jzou1115
# -t 1-1:1
# -t 2-2206:1
#$ -t 2151-2151:1

source ~/.bashrc

ind=$(( SGE_TASK_ID - 1 ))

outdir=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019_Summary_Statistics/permutations
python compute_pval_corrected_perm.py $ind $outdir
