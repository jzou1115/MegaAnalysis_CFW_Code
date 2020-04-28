#!/bin/sh
#$ -N fit
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=5G,highp,h_rt=23:59:59
#$ -m a
#$ -M jzou1115
#$ -t 1-33:1

pheno=`sed -n ${SGE_TASK_ID}p phenotypes.txt`

input=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019_Summary_Statistics/OX/${pheno}.ox.assoc.txt.gz
output=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019_Summary_Statistics/OX/${pheno}.ox.assoc_qqplot.png

Rscript plot_qqplot.R $input $output
