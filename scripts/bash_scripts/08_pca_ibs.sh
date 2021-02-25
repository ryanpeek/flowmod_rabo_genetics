#!/bin/bash

#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J pcaIB
#SBATCH -e slurms/08_pca_ibs.%j.err
#SBATCH -o slurms/08_pca_ibs.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=5-20:00:00

set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

# run script with
#       sbatch --mem MaxMemPerNode 07_pca_ibs.sh bamlistNAME

# set up directories
pop="SOMM446"
#code_dir="/home/rapeek/projects/SEQS/${pop}/code"
data_dir="/home/rapeek/projects/SEQS/${pop}"

cd ${data_dir}

# make results folder
mkdir -p results_pca

# move to bams folder
cd ${data_dir}/bams

echo "Starting Job: "
date

# specify bamlist
bamlist=$1

nInd=$(wc ${data_dir}/bamlists/${bamlist}.bamlist | awk '{print $1}')
#minInd=$[$nInd/2]
minInd=1

# ANGSD Single Read Sampling approach PCA/MDS (http://www.popgen.dk/angsd/index.php/PCA_MDS)

angsd -bam ${data_dir}/bamlists/${bamlist}.bamlist -out ${data_dir}/results_pca/${bamlist}_genos -doIBS 1 -doCounts 1 -doMajorMinor 1 -doMaf 1 -minMaf 0.05 -minInd ${minInd} -minMapQ 30 -minQ 20 -SNP_pval 1e-9 -makeMatrix 1 -doCov 1 -GL 1 -nThreads 16 -sites ${data_dir}/bait_lengths.txt -minFreq 0.05

