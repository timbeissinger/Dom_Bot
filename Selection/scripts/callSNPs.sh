#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/Selection
#SBATCH -J callSNPs
#SBATCH -o slurm-log/callSNPs_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/callSNPs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
set -e
set -u

module load samtools

# This script will call SNPs in the TILs, BKNs, and Tripsicum
# 8-13-2013


### Set values
output=SNPs
ref=/group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa
snpCallList=INS/snpCallList.txt
maizefile=/home/beissing/Dom_Bot_Git/WholeGenome/INS/BKN_list.txt
minMapQ=30
minQ=20

### Run samtools on teosinte
echo "samtools mpileup -u -f $ref -b $snpCallList -q $minMapQ -Q $minQ -C50 -D -S | bcftools view -bvcg - > $output/all.raw.bcf "
samtools mpileup -u -f $ref -b $snpCallList -q $minMapQ -Q $minQ -C50 -D -S | bcftools view -bvcg - > $output/all.raw.bcf 
echo
echo " bcftools view all.raw.bcf | vcfutils.pl varFilter -D100 > $output/all.flt.vcf  "
bcftools view all.raw.bcf | vcfutils.pl varFilter -D100 > $output/all.flt.vcf 


