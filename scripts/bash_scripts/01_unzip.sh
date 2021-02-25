#!/bin/bash
#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J unzip
#SBATCH -e slurms/01_unzip.err
#SBATCH -o slurms/01_unzip.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=4-20:00:00

set -e
set -x

rawseq=$1

# CREATE DIR IN SEQS
mkdir -p /home/rapeek/projects/SEQS/${rawseq}/raw
mkdir -p /home/rapeek/projects/SEQS/${rawseq}/fastq

# SYMLINK RAW DATA
ln -s /home/millermr/UCDavis/${rawseq}/*fq.gz /home/rapeek/projects/SEQS/${rawseq}/raw/

# if already split use this line
#cd /home/rapeek/projects/SEQS/${rawseq}/fastq/
#cp /home/millermr/UCDavis/${rawseq}/*.fastq.gz /home/rapeek/projects/SEQS/${rawseq}/fastq/

#  This script will unzip all raw RAD/RAPTURE data
#       Run with
#               sbatch --mem MaxMemPerNode 01_unzip.sh SOMM504

# UNZIP DATA
cd /home/rapeek/projects/SEQS

date

for i in ${rawseq}
do
cd ${i}/raw
for file in *.gz
do
newname=$(basename $file .gz)
echo "unzipping " $file " to " $newname
gunzip -f $file
#chmod a=r $newname
done
cd ../
done

date

