#!/bin/sh
#$ -N locuszoom
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=15G,highp,h_rt=48:59:59
#$ -m a
#$ -M jzou1115
# -t 1-1:1
#$ -t 15-15:1
source ~/.bashrc

i=$SGE_TASK_ID
#for i in {1..7}; do 
python locuszoom_input_NA.py $i

#done
