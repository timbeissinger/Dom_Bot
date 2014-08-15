#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/Selection
#SBATCH -J callSNPsAngsd
#SBATCH -o slurm-log/callSNPsAngsd_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/callSNPsAngsd_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
set -e
set -u

# This script will call SNPs in the TILs, BKNs, and Tripsicum, using Angsd
# 8-13-2013


### Set values
angsdir=/home/beissing/bin/angsd0.610
output=SNPs
ref=/group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
snpCallList=INS/snpCallList.txt
maizefile=/home/beissing/Dom_Bot_Git/WholeGenome/INS/BKN_list.txt
regionfile="/home/beissing/Dom_Bot_Git/WholeGenome/INS/wholeGenomeRegionFile.txt"
minMapQ=30
minQ=20
glikehood=1
minMapQ=30
cpu=32
SNP_pval=1e-6

### Run angsd on bams to call snps

echo "$angsdir/angsd -bam $snpCallList -GL $glikehood -out $output/angsd_snps -doMaf 1 -doMajorMinor 1 -doGeno 5 -doPost 1 -postCutoff 0.95 -minMapQ $minMapQ -minQ 20 -rf $regionfile -P $cpu -SNP_pval $SNP_pval"
echo 
$angsdir/angsd -bam $snpCallList -GL $glikehood -out $output/angsd_snps -doMaf 1 -doMajorMinor 1 -doGeno 5 -doPost 1 -postCutoff 0.95 -minMapQ $minMapQ -minQ 20 -rf $regionfile -P $cpu -SNP_pval $SNP_pval



