#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -o /home/beissing/Dom_Bot_Git/slurm-log/index_merged_bams.out
#SBATCH -e /home/beissing/Dom_Bot_Git/slurm-log/index_merged_bams.err
#SBATCH -J index_merged_bams
#SBATCH -p bigmem
set -e
set -u

# re-index previously merged bam files

for i in $( ls /group/jrigrp3/bottleneckProject/mergedBams/*.bam ); do
echo $i;
samtools index "$i";
done;
