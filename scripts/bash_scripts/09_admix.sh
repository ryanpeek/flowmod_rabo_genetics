#!/bin/bash

#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J ngsAdmix
#SBATCH -e slurms/09_admix.%j.err
#SBATCH -o slurms/09_admix.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=5-20:00:00

set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

# run script with
#       sbatch --mem MaxMemPerNode 09_ngsAdmix.sh bamlistNAME kGroups

# set up directories
pop="SOMM446"
#code_dir="/home/rapeek/projects/SEQS/${pop}/code"
data_dir="/home/rapeek/projects/SEQS/${pop}"
ref=final_contigs_300.fa # de novo reference alignment

cd ${data_dir}

# make results folder
mkdir -p results_admix

# move to bams folder
cd ${data_dir}/bams

echo "Starting Job: "
date

# specify bamlist
bamlist=$1

# specify k groups
k1=$2

nInd=$(wc ${data_dir}/bamlists/${bamlist}.bamlist | awk '{print $1}')
#minInd=$[$nInd/2]
minInd=3

# PCAngsd (http://www.popgen.dk/software/index.php/PCAngsd)

## First generate genotype likelihoods in Beagle format using angsd
# angsd -bam ${data_dir}/bamlists/${bamlist}.bamlist -out ${data_dir}/results_pca/${bamlist}_genolikes -doGlf 2 -doMajorMinor 1 -minMaf 0.05 -doMaf 2 -minInd ${minInd} -minMapQ 30 -minQ 20 -SNP_pval 1e-9 -GL 1 -nThreads 16 -ref $ref -sites ${data_dir}/bait_lengths.txt

# calculate admixture
x=2
while [ $x -le $k1 ]
do
	~/bin/NGSadmix -likes ${data_dir}/results_pca/${bamlist}_genolikes.beagle.gz -K $x -o ${data_dir}/results_admix/${bamlist}_admix_k${x}
	x=$(( $x + 1 ))
done
echo "All finished!"

