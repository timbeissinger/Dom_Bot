#!/bin/bash -l
#SBATCH -D /home/beissing/Dom_Bot_Git
#SBATCH -J 2dsfs
#SBATCH -o slurm-log/2dsfs.out
#SBATCH -p bigmem
#SBATCH -e slurm-log/2dsfs.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
set -e
set -u

# Timothy M. Beissinger
# 7-3-2014

# first read in the two arguments, which should be the two populations that
# we wish to compute the 2dsfs for.
pop1=$1
pop2=$2

echo $pop1
echo $pop2
echo

# Next, we set initial values
angsdir=/home/beissing/bin/angsd0.609
outputdir=/home/beissing/Dom_Bot_Git/TEMP
nIndPop1=$( wc -l DATA/LISTS/"$pop1"_list.txt | cut -f 1 -d " " )
nIndPop2=$( wc -l DATA/LISTS/"$pop2"_list.txt | cut -f 1 -d " " )
nPop1=$( expr 2 \* $nIndPop1 )
nPop2=$( expr 2 \* $nIndPop2 )
glikehood=1
minMapQ=30
cpu=32
range="10:1-10000000"

# Script will run ANGSD and calculate 2d sfs for pop1 and pop2 hapmap files
# have already run ANGSD for pop1 and pop2 individually, output of interest is
# in "$pop1".saf.pos.gz and "pop2".saf.pos.gz -- we need to know which sites have
# coverage in both populations.
# NOTICE: IF INDIVIDUAL RUNS WERE ON A SUBSET OF POSITIONS, THIS SCRIPT
# SHOULD BE RUN ON THE SAME SUBSET (probably has to be...)

# Next, extract the compressed files

#command1="gunzip "$outputdir"/"$pop1".saf.pos.gz "$outputdir""
#command2="gunzip "$outputdir"/"$pop2".saf.pos.gz "$outputdir""
#echo $command1
echo
#echo $command2
echo
#$command1
#$command2

# Now we find the positions that occur in both populations using the
# uniq POSIX program

command3=" "$outputdir"/"$pop1".saf.pos "$outputdir"/"$pop2".saf.pos|sort|uniq -d >"$outputdir"/intersect."$pop1"."$pop2".txt"
echo cat $command3
echo 
cat "$outputdir"/"$pop1".saf.pos "$outputdir"/"$pop2".saf.pos|sort|uniq -d >"$outputdir"/intersect."$pop1"."$pop2".txt

# Now redo angsd sample allele frequency calculation by conditioning on
# the sites that occur in both populations.

command4=" -bam DATA/LISTS/"$pop1"_list.txt -out "$outputdir"/"$pop1"conditioned -doMajorMinor 1 -doMaf 1 -indF DATA/INBREEDING/"$pop1".indF -doSaf 1 -uniqueOnly 0 -anc DATA/TRIP/TRIP.fa.gz -minMapQ $minMapQ -minQ 20 -nInd $nIndPop1 -baq 1 -ref /home/beissing/GENOMES/Zea_mays.AGPv3.22.dna.genome.fa -GL $glikehood -P $cpu -r $range -sites "$outputdir"/intersect."$pop1"."$pop2".txt"
echo "$angsdir"/angsd $command4
echo 
"$angsdir"/angsd $command4

command5=" -bam DATA/LISTS/"$pop2"_list.txt -out "$outputdir"/"$pop2"conditioned -doMajorMinor 1 -doMaf 1 -indF DATA/INBREEDING/"$pop2".indF -doSaf 1 -uniqueOnly 0 -anc DATA/TRIP/TRIP.fa.gz -minMapQ $minMapQ -minQ 20 -nInd $nIndPop2 -baq 1 -ref /home/beissing/GENOMES/Zea_mays.AGPv3.22.dna.genome.fa -GL $glikehood -P $cpu -r $range -sites "$outputdir"/intersect."$pop1"."$pop2".txt"
echo "$angsdir"/angsd $command5
echo 
"$angsdir"/angsd $command5



# Now we estimate the joint sfs using the realSFS program
command6=" 2dsfs "$outputdir"/"$pop1"conditioned.saf "$outputdir"/"$pop2"conditioned.saf $nPop1 $nPop2 -P $cpu"
echo "$angsdir"/misc/realSFS $command6 " > SFS/2dsfs."$pop1"."$pop2".sfs"
echo 
"$angsdir"/misc/realSFS $command6  > SFS/2dsfs."$pop1"."$pop2".sfs