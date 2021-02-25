#!/bin/bash
#SBATCH -J split_plate
#SBATCH -e slurms/00_split_plate.%j.err
#SBATCH -o slurms/00_split_plate.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=1-20:00:00

#set -e # exits upon failing command
#set -v # verbose -- all lines
#set -x # trace of all commands after expansion before execution

##This script is set grep out the 6bp 96TruSeq plus blackBarcodes - 6bp??+8bp
##run with sbatch 00_split_novoseq.bash raw_fq/SEQ1 raw_fq/SEQ2 SOMM511
mkdir -p slurms


seq1=$1
seq2=$2
somm=$3


# original used both bc but if swapping an issue, use just first 6:
grep --no-group-separator -A3 ':TGCGCT' ${seq1} > ${somm}_TGCGCT_R1.fastq
grep --no-group-separator -A3 ':TGCGCT' ${seq2} > ${somm}_TGCGCT_R2.fastq
grep --no-group-separator -A3 ':GACCCT' ${seq1} > ${somm}_GACCCT_R1.fastq
grep --no-group-separator -A3 ':GACCCT' ${seq2} > ${somm}_GACCCT_R2.fastq
grep --no-group-separator -A3 ':ACACCT' ${seq1} > ${somm}_ACACCT_R1.fastq
grep --no-group-separator -A3 ':ACACCT' ${seq2} > ${somm}_ACACCT_R2.fastq
grep --no-group-separator -A3 ':CCTACT' ${seq1} > ${somm}_CCTACT_R1.fastq
grep --no-group-separator -A3 ':CCTACT' ${seq2} > ${somm}_CCTACT_R2.fastq
grep --no-group-separator -A3 ':TTAACT' ${seq1} > ${somm}_TTAACT_R1.fastq
grep --no-group-separator -A3 ':TTAACT' ${seq2} > ${somm}_TTAACT_R2.fastq
grep --no-group-separator -A3 ':TACTAT' ${seq1} > ${somm}_TACTAT_R1.fastq
grep --no-group-separator -A3 ':TACTAT' ${seq2} > ${somm}_TACTAT_R2.fastq
grep --no-group-separator -A3 ':TTTGAT' ${seq1} > ${somm}_TTTGAT_R1.fastq
grep --no-group-separator -A3 ':TTTGAT' ${seq2} > ${somm}_TTTGAT_R2.fastq
grep --no-group-separator -A3 ':GTGCAT' ${seq1} > ${somm}_GTGCAT_R1.fastq
grep --no-group-separator -A3 ':GTGCAT' ${seq2} > ${somm}_GTGCAT_R2.fastq
grep --no-group-separator -A3 ':CAACAT' ${seq1} > ${somm}_CAACAT_R1.fastq
grep --no-group-separator -A3 ':CAACAT' ${seq2} > ${somm}_CAACAT_R2.fastq

#grep --no-group-separator -A3 '+GGCATACT' ${seq1} > ${somm}_GGCATACT_R1.fastq
#grep --no-group-separator -A3 '+GGCATACT' ${seq2} > ${somm}_GGCATACT_R2.fastq
#grep --no-group-separator -A3 '+CTTCACCA' ${seq1} > ${somm}_CTTCACCA_R1.fastq
#grep --no-group-separator -A3 '+CTTCACCA' ${seq2} > ${somm}_CTTCACCA_R2.fastq
#grep --no-group-separator -A3 '+GAGCAGTA' ${seq1} > ${somm}_GAGCAGTA_R1.fastq
#grep --no-group-separator -A3 '+GAGCAGTA' ${seq2} > ${somm}_GAGCAGTA_R2.fastq
#grep --no-group-separator -A3 '+GTAGGAGT' ${seq1} > ${somm}_GTAGGAGT_R1.fastq
#grep --no-group-separator -A3 '+GTAGGAGT' ${seq2} > ${somm}_GTAGGAGT_R2.fastq
#grep --no-group-separator -A3 '+CTGGAGTA' ${seq1} > ${somm}_CTGGAGTA_R1.fastq
#grep --no-group-separator -A3 '+CTGGAGTA' ${seq2} > ${somm}_CTGGAGTA_R2.fastq
#grep --no-group-separator -A3 '+CTAGGTGA' ${seq1} > ${somm}_CTAGGTGA_R1.fastq
#grep --no-group-separator -A3 '+CTAGGTGA' ${seq2} > ${somm}_CTAGGTGA_R2.fastq
#grep --no-group-separator -A3 '+GCAAGATC' ${seq1} > ${somm}_GCAAGATC_R1.fastq
#grep --no-group-separator -A3 '+GCAAGATC' ${seq2} > ${somm}_GCAAGATC_R2.fastq
#grep --no-group-separator -A3 '+TCAGGCTT' ${seq1} > ${somm}_TCAGGCTT_R1.fastq
#grep --no-group-separator -A3 '+TCAGGCTT' ${seq2} > ${somm}_TCAGGCTT_R2.fastq
#grep --no-group-separator -A3 '+CGTGTACT' ${seq1} > ${somm}_CGTGTACT_R1.fastq
#grep --no-group-separator -A3 '+CGTGTACT' ${seq2} > ${somm}_CGTGTACT_R2.fastq
