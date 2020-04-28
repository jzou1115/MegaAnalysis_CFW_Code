#!/bin/sh
#$ -N kinship
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=10G,highp,h_rt=23:59:59
#$ -m ae
#$ -M jzou1115
#$ -t 15-19:1

source ~/.bashrc
chrm=$SGE_TASK_ID
dosages=/u/home/j/jzou1115/project-zarlab/CFW/MegaAnalysis_CFW_2019/genotypes/oxinfo90/combined.chr${chrm}.oxinfo90.dosages.gz
out=/u/home/j/jzou1115/project-zarlab/CFW/confidence_intervals/output/grms/combined.chr${chrm}.oxinfo90.kinship.RData

Rscript make_kinship_matrices.R $dosages $out

