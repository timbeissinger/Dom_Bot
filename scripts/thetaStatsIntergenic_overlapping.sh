#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J thetaStatsIntergenic_overlapping
#SBATCH -o slurm-log/thetaStatsIntergenic_overlapping_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/thetaStatsIntergenic_overlapping_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# script to compute thetas and Tajima's D for already-run
# sfs.

### Set variables
angsdir=/home/beissing/bin/angsd0.609

taxon=$1
append=$2
nInd=$( wc -l DATA/LISTS/"$taxon"_list.txt | cut -f 1 -d " " )
n=$( expr 2 \* $nInd )
minperc=0.8
minInd=$( printf "%.0f" $(echo "scale=2;$minperc*$nInd" | bc))
glikehood=1
minMapQ=30
cpu=32
regionfile="/home/beissing/Dom_Bot_Git/intergenic/intergenicRegionFile_chr10.txt"
windowsize=1000
stepsize=500

### Echo params
echo "Taxon = "$taxon""
echo "File append = "$append""

### Calculate sfs on only overlapping sites
command0="/home/beissing/Dom_Bot_Git/TEMP/"$taxon"_intergenic_10_conditioned.saf $n -P $cpu"
echo $command0
$angsdir/misc/realSFS $command0 > SFS/"$taxon"_intergenic_10_conditioned.sfs

### Calculate thetas for each site
command1="-bam DATA/LISTS/"$taxon"_list.txt -out THETAS/"$taxon"_"$append" -pest SFS/"$taxon"_intergenic_10_conditioned.sfs -indF DATA/INBREEDING/"$taxon".indF -doSaf 2 -uniqueOnly 0 -anc DATA/TRIP/TRIP.fa.gz -minMapQ $minMapQ -minQ 20 -nInd $nInd -minInd $minInd -baq 1 -ref /group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa -GL $glikehood -P $cpu -rf $regionfile -doThetas 1 -doMajorMinor 1 -doMaf 1 -sites /home/beissing/Dom_Bot_Git/TEMP/intersect.TIL.BKN_intergenic_10.txt"
echo $command1
$angsdir/angsd $command1

### Estimate Tajima's D and other statistics
command2="make_bed /home/beissing/Dom_Bot_Git/THETAS/"$taxon"_"$append".thetas.gz "
echo $command2
$angsdir/misc/thetaStat $command2

command3=" do_stat  /home/beissing/Dom_Bot_Git/THETAS/"$taxon"_"$append".thetas.gz -nChr $n -win $windowsize -step $stepsize "
echo $command3
$angsdir/misc/thetaStat $command3