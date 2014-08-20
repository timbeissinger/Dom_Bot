#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/Selection
#SBATCH -J makeFasta
#SBATCH -o /home/beissing/Dom_Bot_Git/slurm-log/makeFasta_%j.out
#SBATCH -p bigmem
#SBATCH -e /home/beissing/Dom_Bot_Git/slurm-log/makeFasta_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
set -e
set -u


# This script will create fasta files for each of the merged bam files, to
# be used later by polydNdS
#8-20-2014

### Set values
angsdir=/home/beissing/bin/angsd0.610
ref=/group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
taxaList=INS/snpCallList.txt

### Run
echo "$angsdir/angsd -bam $taxaList -doFasta 1 -out /group/jrigrp3/bottleneckProject/fasta"

$angsdir/angsd -bam $taxaList -doFasta 1 -out /group/jrigrp3/bottleneckProject/fasta
