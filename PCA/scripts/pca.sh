#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J pca
#SBATCH -o slurm-log/pca_%j.out
#SBATCH -p serial
#SBATCH -e slurm-log/pca_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

# Script to compute PCA from angsd output
taxon=BKN
genosfile=BKN_begin_10
ngsdir=/home/beissing/bin/ngsTools/ngsPopGen
datadir=/home/beissing/TEMP


nInd=$( wc -l DATA/LISTS/"$taxon"_list.txt | cut -f 1 -d " " )



# Compute PCA
command1="-probfile "$datadir"/"$genosfile".geno.gz -outfile /home/beissing/Dom_Bot_Git/PCA/"$taxon".covar1 -nind $nInd -nsites 10000 -call 0  -minmaf 0.05"
echo $command1
$ngsdir/ngsCovar $command1