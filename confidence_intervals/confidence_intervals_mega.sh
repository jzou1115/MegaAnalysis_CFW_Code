#!/bin/sh
#$ -N ci
#$ -cwd
#$ -o output
#$ -e error
#$ -V
# -l h_data=20G,exclusive,highp,h_rt=96:59:59
#$ -l h_data=20G,exclusive,highp,h_rt=96:59:59
# -l h_data=30G,highp,h_rt=96:59:59
# -pe shared 5
#$ -pe dc\* 5
#$ -m abe
#$ -M jzou1115

source ~/.bashrc

Rscript confidence_intervals_mega.R
