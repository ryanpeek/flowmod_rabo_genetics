#!/bin/bash
#SBATCH -J align
#SBATCH -e slurms/04_align_%j.err
#SBATCH -o slurms/04_align_%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=4-20:00:00

#set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

#  This script will align all raw fasta files. Run with:
#       sbatch 04_align.sh somm barcode
#  assume running script while located in fastq dir

# INPUTS
somm=$1 # the full path to the fastqs 

# switch to dir
echo "changing dir"
cd ${somm}

#if [ ! -d "/${barcode}" ] && echo "Directory ${barcode} DOES NOT exist..."; then
#        cd ${barcode}
#        fi

# make a bamlist for alignment: 3 columns, RA, RB output_name (something without AB)
ls *RA* > list_A; ls *RB* > list_B; ls *RA* | sed "s/\_RA//g" | sed "s/\.fastq//g" > list_noAB; paste list_* > align_list

rm list* # remove temp lists

# set list
list="align_list"

# set reference alignment
#ref=$2 # give reference alignment
ref="/home/rapeek/projects/SEQS/final_contigs_300.fa"

wc=$(wc -l ${list} | awk '{print $1}')

x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1, $2, $3}')   
	set -- $var
	c1=$1
	c2=$2
	c3=$3

	echo "#!/bin/bash -l
#SBATCH -o slurms/04_bwa_%j.out
#SBATCH -p high

bwa mem $ref ${c1} ${c2} | samtools view -Sb - | samtools sort - -o ${c3}.sort.bam
samtools view -f 0x2 -b ${c3}.sort.bam | samtools rmdup - ${c3}.sort.flt.bam" > ${c3}.sh

	sbatch --mem=8G -t 8:00:00 ${c3}.sh

	x=$(( $x + 1 ))

done

