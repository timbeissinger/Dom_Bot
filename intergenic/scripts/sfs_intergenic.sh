#!/bin/bash -l

#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J sfs_intergenic
#SBATCH -o slurm-log/sfs_intergenic_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/sfs_intergenic_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

# script to run ANGSD on hapmap2 bam files
module load angsd



angsdir=/home/beissing/bin/angsd0.609

#windowsize=1000
#step=500
taxon=$1
nInd=$( wc -l DATA/LISTS/"$taxon"_list.txt | cut -f 1 -d " " )
n=$( expr 2 \* $nInd )
minperc=0.8
minInd=$( printf "%.0f" $(echo "scale=2;$minperc*$nInd" | bc))
glikehood=1
minMapQ=30
cpu=32
regionfile="/home/beissing/Dom_Bot_Git/intergenic/intergenicRegionFile_chr10.txt"
#range="10:1-10000000"

# (estimate an SFS)
# -bam list of paths to bamfiles you want to use
# -out output file (prior for SFS I believe)
# -doSaf estimates SFS, do 2 if you have inbreeding coefficients (-indF)
# -uniqueOnly 1=use only uniquely mapping reads
# -anc ancestral sequence (see trip.sh script for how this is generated)
# -minMapQ 30 minimum mapping quality of reads to accepy
# -minQ 20 minimum bp quality
# -setMaxDepth 20 sets max depth to accept -- useful to deal with highly repetitive regions
# -baq 1=realign locally (I think)
# -GL $glikehood 1 is samtools, 2 is GATK, 3 SOAPsnp 4 SYK
# -rf path to intergenic text file
# -P 8 use 8 threads
# -indF individiual inbreeding coefficient. for inbred lines just make a files of "1" on each line for each bamfile. otherwise use ngsF to estimate (see inbreeding.sh script)

command1="-bam DATA/LISTS/"$taxon"_list.txt -out /home/beissing/Dom_Bot_Git/TEMP/"$taxon"_intergenic_10 -doMajorMinor 1 -doMaf 1 -indF DATA/INBREEDING/"$taxon".indF -doSaf 2 -uniqueOnly 0 -anc DATA/TRIP/TRIP.fa.gz -minMapQ $minMapQ -minQ 20 -nInd $nInd -minInd $minInd -baq 1 -ref /group/jrigrp3/bottleneckProject/genomes/Zea_mays.AGPv3.22.dna.genome.fa -GL $glikehood -P $cpu -rf $regionfile"
echo $command1
$angsdir/angsd $command1


# not clear to me how to run folded, as -fold option seems to be deprecated?
# temp/"$taxon"_pest.saf output file from above run; prior on SFS?
# $n number of chromosomes; 2 x number of inds for diploids
# results/"$taxon"_pest.em.ml This output is the final estimated SFS
	# the file will be nat. log probabilities of the value of the SFS from 0:n
	# so if n=10, there will be 11 numbers.  To plot the SFS for polymorphic sites only, ignore the first and last numbers. e.g. for teosinte I get:
	# -0.133730 -3.724029 -4.246469 -4.981319 -5.453217 -5.803669 -6.076224 -6.330416 -6.501992 -6.713127 -6.882129 -6.970549 -7.289374 -7.434923 -7.308903 -7.057695 -7.457825 -7.740251 -7.665521 -7.683324 -7.788163 -7.702094 -7.562837 -7.491339 -7.416449 -7.364919 -7.107873 -6.870063 -6.458559 -6.044445 -2.994086
	# which corresponds to exp(-0.13)~0.9 or 90% of sites are fixed for ancestral allele, and exp(-2.994086) or ~5% are fixed for derived allele. 
	# remaining 5% are polymorphic
command2="/home/beissing/Dom_Bot_Git/TEMP/"$taxon"_intergenic_10.saf $n -P $cpu" 
echo $command2
$angsdir/misc/realSFS $command2 > SFS/"$taxon"_intergenic_10.sfs
