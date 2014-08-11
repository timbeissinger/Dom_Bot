#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J indexRefs
#SBATCH -o slurm-log/indexRefs_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/indexRefs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
set -e
set -u

#echo "step 1"
#samtools faidx /group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
echo
echo "step 2"
samtools faidx /group/jrigrp3/bottleneckProject/genomes/TRIP.fa