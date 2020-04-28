#!/bin/sh
#$ -N ci3
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=20G,highp,h_rt=48:59:59
#$ -m abe
#$ -M jzou1115

source ~/.bashrc

Rscript confidence_intervals_mega3.R
