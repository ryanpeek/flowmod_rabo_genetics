#!/bin/bash
#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J compressFQ
#SBATCH -e 18_compress_fastq.%j.err
#SBATCH -o 18_compress_fastq.%j.out
#SBATCH -c 20
#SBATCH -p med
#SBATCH --time=1-20:00:00

set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

dir=$1 # give a directory to compress

# /home/rapeek/projects/SEQS/SOMM446/fastq

cd $1
echo Compressing SOMM R* fastq files on $(date | awk '{print $4 " " $3 $2 $6}')
lzma *.fastq
#gzip *.fastq
# -v to see compression, default is 6
# -r to compress all files in a dir

echo Compression complete at $(date | awk '{print $4 " " $3 $2 $6}')
