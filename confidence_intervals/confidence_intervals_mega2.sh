#!/bin/sh
#$ -N ci2
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=4G,highp,h_rt=48:59:59
#$ -pe shared 10
#$ -m abe
#$ -M jzou1115

source ~/.bashrc

Rscript confidence_intervals_mega2.R
