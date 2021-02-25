#!/bin/bash
#SBATCH -J del_fq
#SBATCH -c 10
#SBATCH -p high
#SBATCH --time=1-20:00:00

set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

# You may need modifications to run it alone.
# run with sbatch -t 1200 -p high 03_split_wells.sh SOMM504 barcode.txt 

# assumes in SOMMXXX/fastq directory where fastq R1/R2 are
fqdir=$1

# list of barcodes of interest
bc1=$2 # ATGGCT, ACACCT, ATGGCT

# now a while loop to read each line of bc1
while IFS= read -r barcode
do
	echo $barcode
	if [ ! -d "/${barcode}" ] && echo "Directory ${barcode} DOES NOT exist..."; then	
	rm ${barcode}/${fqdir}_${barcode}_R1.fastq
	rm ${barcode}/${fqdir}_${barcode}_R2.fastq
	echo "done"
	fi

done < ${bc1}

