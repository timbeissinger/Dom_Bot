#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -o /home/beissing/Dom_Bot_Git/slurm-log/bam_merge_out.txt
#SBATCH -e /home/beissing/Dom_Bot_Git/slurm-log/bam_merge_err.txt
#SBATCH -J bam_merge
#SBATCH -p bigmem
set -e
set -u

# merge bam files to single individuals 

for i in $( ls /home/beissing/BAMLINKS/ | cut -f 1 -d "_" | sort -n | uniq ); do

#if not already merged
if [ ! -f /home/beissing/BAMLINKS/"$i"_merged.bam ]
then

	if [ $( ls -1 /home/beissing/BAMLINKS/$i* | wc -l ) -gt 1 ]
	then
		samtools merge -r /home/beissing/BAMLINKS/"$i"_merged.bam $( ls /home/beissing/BAMLINKS/$i* | perl -ne '{BEGIN} chomp; print "$_\t"; while(<>){ chomp; print "$_\t";}; {END} chomp; print $_;')
		samtools index /home/beissing/BAMLINKS/"$i"_merged.bam
	else
		mv $( ls -1 /home/beissing/BAMLINKS/$i*  ) /home/beissing/v3_bams_bwamem/"$i"_merged.bam;
		samtools index /home/beissing/BAMLINKS/"$i"_merged.bam
	fi;
else 
 	if [ ! -f /home/beissing/BAMLINKS/"$i"_merged.bam.bai ]
	then
		samtools index /home/beissing/BAMLINKS/"$i"_merged.bam
	fi;
fi;

done;
