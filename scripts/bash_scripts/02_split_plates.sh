#!/bin/bash
#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J split
#SBATCH -e slurms/02_split.%j.err
#SBATCH -o slurms/02_split.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=1-20:00:00

#set -e # exits upon failing command
#set -v # verbose -- all lines
#set -x # trace of all commands after expansion before execution

mkdir -p slurms

seq1=$1
seq2=$2
somm=$3

grep --no-group-separator -A 3 ":TCAGTT" ${seq1} > ${somm}_TCAGTT_R1.fastq
grep --no-group-separator -A 3 ":TCAGTT" ${seq2} > ${somm}_TCAGTT_R2.fastq

grep --no-group-separator -A 3 ":CGTCTT" ${seq1} > ${somm}_CGTCTT_R1.fastq
grep --no-group-separator -A 3 ":CGTCTT" ${seq2} > ${somm}_CGTCTT_R2.fastq

grep --no-group-separator -A 3 ":CTCATT" ${seq1} > ${somm}_CTCATT_R1.fastq
grep --no-group-separator -A 3 ":CTCATT" ${seq2} > ${somm}_CTCATT_R2.fastq

grep --no-group-separator -A 3 ":TCTTGT" ${seq1} > ${somm}_TCTTGT_R1.fastq
grep --no-group-separator -A 3 ":TCTTGT" ${seq2} > ${somm}_TCTTGT_R2.fastq

grep --no-group-separator -A 3 ":CGCTGT" ${seq1} > ${somm}_CGCTGT_R1.fastq
grep --no-group-separator -A 3 ":CGCTGT" ${seq2} > ${somm}_CGCTGT_R2.fastq

grep --no-group-separator -A 3 ":GTATGT" ${seq1} > ${somm}_GTATGT_R1.fastq
grep --no-group-separator -A 3 ":GTATGT" ${seq2} > ${somm}_GTATGT_R2.fastq

grep --no-group-separator -A 3 ":CATGGT" ${seq1} > ${somm}_CATGGT_R1.fastq
grep --no-group-separator -A 3 ":CATGGT" ${seq2} > ${somm}_CATGGT_R2.fastq

grep --no-group-separator -A 3 ":GGGGGT" ${seq1} > ${somm}_GGGGGT_R1.fastq
grep --no-group-separator -A 3 ":GGGGGT" ${seq2} > ${somm}_GGGGGT_R2.fastq

grep --no-group-separator -A 3 ":AAGAGT" ${seq1} > ${somm}_AAGAGT_R1.fastq
grep --no-group-separator -A 3 ":AAGAGT" ${seq2} > ${somm}_AAGAGT_R2.fastq

grep --no-group-separator -A 3 ":GCCAGT" ${seq1} > ${somm}_GCCAGT_R1.fastq
grep --no-group-separator -A 3 ":GCCAGT" ${seq2} > ${somm}_GCCAGT_R2.fastq

grep --no-group-separator -A 3 ":ATGGCT" ${seq1} > ${somm}_ATGGCT_R1.fastq
grep --no-group-separator -A 3 ":ATGGCT" ${seq2} > ${somm}_ATGGCT_R2.fastq

grep --no-group-separator -A 3 ":TGCGCT" ${seq1} > ${somm}_TGCGCT_R1.fastq
grep --no-group-separator -A 3 ":TGCGCT" ${seq2} > ${somm}_TGCGCT_R2.fastq

grep --no-group-separator -A 3 ":GACCCT" ${seq1} > ${somm}_GACCCT_R1.fastq
grep --no-group-separator -A 3 ":GACCCT" ${seq2} > ${somm}_GACCCT_R2.fastq

grep --no-group-separator -A 3 ":ACACCT" ${seq1} > ${somm}_ACACCT_R1.fastq
grep --no-group-separator -A 3 ":ACACCT" ${seq2} > ${somm}_ACACCT_R2.fastq

grep --no-group-separator -A 3 ":CCTACT" ${seq1} > ${somm}_CCTACT_R1.fastq
grep --no-group-separator -A 3 ":CCTACT" ${seq2} > ${somm}_CCTACT_R2.fastq

grep --no-group-separator -A 3 ":CGATTG" ${seq1} > ${somm}_CGATTG_R1.fastq
grep --no-group-separator -A 3 ":CGATTG" ${seq2} > ${somm}_CGATTG_R2.fastq

