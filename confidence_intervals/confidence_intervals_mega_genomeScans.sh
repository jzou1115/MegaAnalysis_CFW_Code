#!/bin/sh
#$ -N ci
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=50G,highp,h_rt=23:59:59
#$ -m abe
#$ -M jzou1115
# -t 3-34:1

source ~/.bashrc

mmf=$1
#mmf=`sed -n ${SGE_TASK_ID}p mm_files.txt`
Rscript confidence_intervals_mega_genomeScans.R $mmf
