#!/bin/bash

module add r/4.0.1

n=$1 # Number of clusters per each simulation
s=$2 # Number of sample splitting repetition

d=data/n${n}_s${s}

echo $d
mkdir -p $d
mkdir -p $d/log
mkdir -p $d/Rdata
cd $d
cp ../../deltas.rds .

### Estimator Simulation ###
jid_estimator=$( sbatch --array=1-50 ../../estimator.sh $n $s | awk '{print $4}' )

### Estimand Simulation ###
jid_estimand=$( sbatch --array=1-100 ../../estimand.sh | awk '{print $4}' )

### Read Simulation result ###
# jid_readresult=$( sbatch --dependency=afterany:$jid_estimator -o simulation.log --wrap="Rscript ../../readresult.R" | awk '{print $4}' )
jid_readresult=$( sbatch --dependency=afterany:$jid_estimator:$jid_estimand -o simulation.log --wrap="Rscript ../../readresult.R" | awk '{print $4}' )

### Clear intermediate files ###
# sbatch --dependency=afterok:$jid_readresult --wrap="rm -r ./Rdata"