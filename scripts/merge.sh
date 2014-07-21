#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -o /home/beissing/Dom_Bot_Git/slurm-log/bam_merge_out.txt
#SBATCH -e /home/beissing/Dom_Bot_Git/slurm-log/bam_merge_err.txt
#SBATCH -J bam_merge
#SBATCH -p bigmem
set -e
set -u

# merge bam files to single individuals 

# set directories
input="/group/jrigrp3/bottleneckProject/v3_bams_bwamem"
output="/group/jrigrp3/bottleneckProject/mergedBams"

for i in $( ls "$input"/ | cut -f 1 -d "_" | sort -n | uniq ); do

#if not already merged
if [ ! -f "$output"/"$i"_merged.bam ]
then

	if [ $( ls -1 "$input"/$i* | wc -l ) -gt 1 ]
	then
		samtools merge -r "$input"/"$i"_merged.bam $( ls "$input"/$i* | perl -ne '{BEGIN} chomp; print "$_\t"; while(<>){ chomp; print "$_\t";}; {END} chomp; print $_;')
		samtools index  "$input"/"$i"_merged.bam
	else
		mv $( ls -1 "$input"/$i*  ) "$input"/"$i"_merged.bam;
		samtools index "$input"/"$i"_merged.bam
	fi;
else 
 	if [ ! -f "$input"/"$i"_merged.bam.bai ]
	then
		samtools index "$input"/"$i"_merged.bam
	fi;
fi;

done;
