#!/bin/sh
#$ -N gemma
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=20G,highp,h_rt=23:59:59
#$ -m a
#$ -M jzou1115
#$ -t 5-37:1

n=$SGE_TASK_ID
python format_gemma_output.py $n
