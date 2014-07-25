#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J thetaStats
#SBATCH -o slurm-log/thetaStats_%j.out
#SBATCH -p serial
#SBATCH -e slurm-log/thetaStats_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# script to compute thetas and Tajima's D for already-run
# sfs.

### Set variables
angsdir=/home/beissing/bin/angsd0.609

taxon=$1
sfs=$2
append=$3
nInd=$( wc -l DATA/LISTS/"$taxon"_list.txt | cut -f 1 -d " " )
n=$( expr 2 \* $nInd )
minperc=0.8
minInd=$( printf "%.0f" $(echo "scale=2;$minperc*$nInd" | bc))
glikehood=1
minMapQ=30
cpu=32
regionfile="/home/beissing/Dom_Bot_Git/genic/genicRegionFile_chr10.txt"

### Calculate thetas for each site
command1="-bam DATA/LISTS/"$taxon"_list.txt -out THETAS/"$taxon"_"$append" -pest SFS/"$sfs" -indF DATA/INBREEDING/"$taxon".indF -doSaf 2 -uniqueOnly 0 -anc DATA/TRIP/TRIP.fa.gz -minMapQ $minMapQ -minQ 20 -nInd $nInd -minInd $minInd -baq 1 -ref /group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa -GL $glikehood -P $cpu -rf $regionfile -doThetas 1 -doMajorMinor 1 -doMaf 1 "
echo $command1
$angsdir/angsd $command1

### Estimate Tajima's D and other statistics
command2="make_bed /home/beissing/Dom_Bot_Git/THETAS/"$taxon"_"$append".thetas.gz "
echo $command2
$angsdir/misc/thetaStat $command2

command3=" do_stat  /home/beissing/Dom_Bot_Git/THETAS/"$taxon"_"$append".thetas.gz -nChr $n  "
echo $command3
$angsdir/misc/thetaStat $command3