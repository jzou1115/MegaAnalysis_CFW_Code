#!/bin/sh
#$ -N preprocess
#$ -cwd
#$ -o output
#$ -e error
#$ -V
#$ -l h_data=20G,highp,h_rt=23:59:59
#$ -m ae
#$ -M jzou1115

python preprocess_gemma_input_files.py
