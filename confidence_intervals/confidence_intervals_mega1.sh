#!/bin/sh
#$ -N ci20
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=40G,highp,h_rt=23:59:59
#$ -m abe
#$ -M jzou1115

source ~/.bashrc
phenotype=$1
mm=/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output2/emma.${phenotype}/emma.${phenotype}.RData
Rscript confidence_intervals_mega1.R $mm
