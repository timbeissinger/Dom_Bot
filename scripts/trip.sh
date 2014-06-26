#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J trip_fasta
#SBATCH -o /home/beissing/Dom_Bot_Git/slurm-log/trip_fasta_out.txt
#SBATCH -p bigmem
#SBATCH -e /home/beissing/Dom_Bot_Git/slurm-log/trip_fasta_err.txt

module load angsd

angsd -i /home/beissing/BAMLINKS/TDD39103_merged.bam -doFasta 1 -out /home/beissing/Dom_Bot_Git/DATA/TRIP
