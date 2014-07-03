#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -o /home/beissing/Dom_Bot_Git/slurm-log/index_reference.out
#SBATCH -e /home/beissing/Dom_Bot_Git/slurm-log/inex_reference.err
#SBATCH -J index_reference
#SBATCH -p bigmem

samtools faidx /home/beissing/GENOMES/Zea_mays.AGPv3.22.dna.genome.fa