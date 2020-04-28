#!/bin/sh
#$ -N locuszoom
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=15G,highp,h_rt=01:59:59
#$ -m a
#$ -M jzou1115
#$ -t 2-39:1

source ~/.bashrc

python wellcome_ci_input.py  $SGE_TASK_ID

