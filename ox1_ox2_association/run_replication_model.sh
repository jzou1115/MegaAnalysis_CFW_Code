#!/bin/sh
#$ -N fit
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=5G,highp,h_rt=23:59:59
#$ -m a
#$ -M jzou1115
# -t 4-4:1
#$ -t 4-37:1


pheno=$SGE_TASK_ID

#for pheno in {4..37}; do 
Rscript create_model_input.R $pheno
#done
