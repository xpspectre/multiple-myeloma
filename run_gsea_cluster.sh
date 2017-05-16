#!/usr/bin/env bash
#PBS -j oe
#PBS -q default
#PBS -r n
#PBS -l nodes=1:ppn=1
#PBS -l mem=8192Mb,vmem=8192Mb
#PBS -l walltime=0:1:00:00
#PBS -N multiple_myeloma_gsea
#PBS -t 1-676

cd /data/kshi/multiple_myeloma/processed/gsea2
./run_${PBS_ARRAYID}.sh
