#!/bin/bash
#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J compressFQ
#SBATCH -e decompress_fastq.%j.err
#SBATCH -o decompress_fastq.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=1-20:00:00

set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

dir=$1 # give a directory to compress

# /home/rapeek/projects/SEQS/SOMM446/fastq

cd $1
echo Decompressing SOMM R* fastq files on $(date | awk '{print $4 " " $3 $2 $6}')

unlzma *.fastq.lzma
#unzip *.fastq.gz
echo Decompression complete at $(date | awk '{print $4 " " $3 $2 $6}')

#echo 'compress as gz'
#gzip -7 *.fastq

#echo "compression complete at $(date | awk '{print $4 " " $3 $2 $6}')

