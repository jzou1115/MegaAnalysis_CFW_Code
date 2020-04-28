#!/bin/sh
#$ -N fit
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=10G,highp,h_rt=23:59:59
#$ -m a
#$ -M jzou1115
# -t 1-1:1
#$ -t 2-33:1


pheno=`sed -n ${SGE_TASK_ID}p phenotypes.txt`

#python preprocess_summary_statistics2.py $pheno



threshold=.000001
summary=../../ox1_ox2_association/model_in/ox1_ox2_${pheno}_90.txt
out=../../ox1_ox2_association/model_out/ox1_ox2_${pheno}_out_${threshold}.txt
Rscript predict_replication2.R $summary $out $threshold --none
