#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git/Selection
#SBATCH -J VEP
#SBATCH -o slurm-log/VEP_%j.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/VEP_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
set -e
set -u

### Identify variant effects for maizeTeo/trip subs
vepDir=/home/beissing/bin/ensembl-tools-release-76/scripts/variant_effect_predictor

echo $HOME

echo "perl $vepDir/variant_effect_predictor.pl --cache --cache_version 23 --genomes --species zea_mays --assembly AGPv3 -i SNPs/TvMT.vep -o SNPs/TvMTeffects.txt"

perl $vepDir/variant_effect_predictor.pl --cache --cache_version 23 --genomes --species zea_mays --assembly AGPv3 -i SNPs/TvMT.vep -o SNPs/TvMTeffects.txt