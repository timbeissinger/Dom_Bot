#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/Selection
#SBATCH -J findSubstitutions
#SBATCH -o slurm-log/findSubstitutions_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/findSubstitutions_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
set -e
set -u

module load R

# This script will run an R script to identify substitutions

R --no-save < scripts/findSubstitutions.R 
