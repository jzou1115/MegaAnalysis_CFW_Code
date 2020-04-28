#!/bin/sh
#$ -N locuszoom
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=10G,highp,h_rt=48:59:59
#$ -m abe
#$ -M jzou1115

source ~/.bashrc

for i in {1..32}; do 
python locuszoom_input.py $i

done
