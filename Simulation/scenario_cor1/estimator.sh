#!/bin/bash
#SBATCH -t 40:00:00
#SBATCH --output=log/estimator_rep-%a.out

n=$1 # Number of clusters per simulation
s=$2 # Number of sample splitting repetition

Rscript "../../estimator.R" \
  -n $n \
  -K 2 \
  -r 100 \
  -s $s