#!/bin/bash
#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J splitSEQ
#SBATCH -e slurms/02_split_RAD.%j.err
#SBATCH -o slurms/02_split_RAD.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=2-20:00:00

set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution

mkdir -p slurms

ddir=$1
# run script with
#	sbatch --mem MaxMemPerNode 02_split_RAD.sh

	#################
	###  Set Up   ###
	#################
	# Change barcodes
	#################

# set up directories
code_dir="/home/rapeek/projects/SEQS/code"
data_dir="/home/rapeek/projects/SEQS/${ddir}/"

cd $data_dir

for i in ${data_dir} #directory
do
barcode="GCCAGT,ATGGCT,TGCGCT,GACCCT" #barcodes
cd fastq

###################
### Split Lanes ###
###################

## split lanes based on index read
$code_dir/BarcodeSplitList3Files.pl ${i}_*_R1_001.fastq ${i}_*_R2_001.fastq ${i}_*_R4_001.fastq $barcode $i

chmod a=r *.fastq

###################
### Split wells ###
###################

# get list of index bcs
ls ${i}_R3_??????.fastq | sed 's/.*R3_//g' | sed 's/.fastq//g' > indexid.list
head indexid.list

# make slurm dir for err / out files
mkdir -p slurms

# split each lane by wells via inline bc
wc=$(wc -l indexid.list | awk '{print $1}') # number of lines
echo "There are" $wc "barcodes to split for" $i
x=1

while [ $x -le $wc ]
do
string="sed -n ${x}p indexid.list"
str=$($string)
var=$(echo $str | awk -Ft '{print $1}')
set -- $var
c1=$1
echo "creating script to split by" $c1
echo "#!/bin/bash" > 04_split_wells.${i}.${c1}.sh
echo "#SBATCH --mail-user=rapeek@ucdavis.edu" >> 04_split_wells.${i}.${c1}.sh
echo "#SBATCH --mail-type=ALL" >> 04_split_wells.${i}.${c1}.sh 
echo "#SBATCH -J ${c1}.${i}" >> 04_split_wells.${i}.${c1}.sh
echo "#SBATCH -e slurms/04_split_wells.${i}.${c1}.%j.err" >> 04_split_wells.${i}.${c1}.sh
echo "#SBATCH -o slurms/04_split_wells.${i}.${c1}.%j.out" >> 04_split_wells.${i}.${c1}.sh
echo "#SBATCH -c 20" >> 04_split_wells.${i}.${c1}.sh
echo "#SBATCH -p high" >> 04_split_wells.${i}.${c1}.sh
echo "#SBATCH --time=1-20:00:00" >> 04_split_wells.${i}.${c1}.sh
echo "" >> 04_split_wells.${i}.${c1}.sh
echo "set -e # exits upon failing command" >> 04_split_wells.${i}.${c1}.sh
echo "set -v # verbose -- all lines" >> 04_split_wells.${i}.${c1}.sh
echo "set -x # trace of all commands after expansion before execution" >> 04_split_wells.${i}.${c1}.sh
echo "" >> 04_split_wells.${i}.${c1}.sh
echo "# This script is run from ${code_dir}/03_split_lanes.sh" >> 04_split_wells.${i}.${c1}.sh
echo "# You may need modifications to run it alone." >> 04_split_wells.${i}.${c1}.sh
echo "" >> 04_split_wells.${i}.${c1}.sh
echo "cd ${data_dir}/fastq" >> 04_split_wells.${i}.${c1}.sh
echo "" >> 04_split_wells.${i}.${c1}.sh
echo "#for file in SOMM44?_*R1_??????.fastq" >> 04_split_wells.${i}.${c1}.sh
echo "#do" >> 04_split_wells.${i}.${c1}.sh
echo "##  put all the same lane in a directory" >> 04_split_wells.${i}.${c1}.sh
echo "mkdir -p ${c1}" >> 04_split_wells.${i}.${c1}.sh
echo "mv *_${c1}.fastq ${c1}/." >> 04_split_wells.${i}.${c1}.sh
echo "cd ${c1}" >> 04_split_wells.${i}.${c1}.sh
echo "" >> 04_split_wells.${i}.${c1}.sh
echo "$code_dir/BCsplitBestRadPE2.pl ${i}_R1_${c1}.fastq ${i}_R3_${c1}.fastq GGACAAGCTATGCAGG,GGAAACATCGTGCAGG,GGACATTGGCTGCAGG,GGACCACTGTTGCAGG,GGAACGTGATTGCAGG,GGCGCTGATCTGCAGG,GGCAGATCTGTGCAGG,GGATGCCTAATGCAGG,GGAACGAACGTGCAGG,GGAGTACAAGTGCAGG,GGCATCAAGTTGCAGG,GGAGTGGTCATGCAGG,GGAACAACCATGCAGG,GGAACCGAGATGCAGG,GGAACGCTTATGCAGG,GGAAGACGGATGCAGG,GGAAGGTACATGCAGG,GGACACAGAATGCAGG,GGACAGCAGATGCAGG,GGACCTCCAATGCAGG,GGACGCTCGATGCAGG,GGACGTATCATGCAGG,GGACTATGCATGCAGG,GGAGAGTCAATGCAGG,GGAGATCGCATGCAGG,GGAGCAGGAATGCAGG,GGAGTCACTATGCAGG,GGATCCTGTATGCAGG,GGATTGAGGATGCAGG,GGCAACCACATGCAGG,GGCAAGACTATGCAGG,GGCAATGGAATGCAGG,GGCACTTCGATGCAGG,GGCAGCGTTATGCAGG,GGCATACCAATGCAGG,GGCCAGTTCATGCAGG,GGCCGAAGTATGCAGG,GGCCGTGAGATGCAGG,GGCCTCCTGATGCAGG,GGCGAACTTATGCAGG,GGCGACTGGATGCAGG,GGCGCATACATGCAGG,GGCTCAATGATGCAGG,GGCTGAGCCATGCAGG,GGCTGGCATATGCAGG,GGGAATCTGATGCAGG,GGGACTAGTATGCAGG,GGGAGCTGAATGCAGG,GGGATAGACATGCAGG,GGGCCACATATGCAGG,GGGCGAGTAATGCAGG,GGGCTAACGATGCAGG,GGGCTCGGTATGCAGG,GGGGAGAACATGCAGG,GGGGTGCGAATGCAGG,GGGTACGCAATGCAGG,GGGTCGTAGATGCAGG,GGGTCTGTCATGCAGG,GGGTGTTCTATGCAGG,GGTAGGATGATGCAGG,GGTATCAGCATGCAGG,GGTCCGTCTATGCAGG,GGTCTTCACATGCAGG,GGTGAAGAGATGCAGG,GGTGGAACAATGCAGG,GGTGGCTTCATGCAGG,GGTGGTGGTATGCAGG,GGTTCACGCATGCAGG,GGACACGAGATGCAGG,GGAAGAGATCTGCAGG,GGAAGGACACTGCAGG,GGAATCCGTCTGCAGG,GGAATGTTGCTGCAGG,GGACACTGACTGCAGG,GGACAGATTCTGCAGG,GGAGATGTACTGCAGG,GGAGCACCTCTGCAGG,GGAGCCATGCTGCAGG,GGAGGCTAACTGCAGG,GGATAGCGACTGCAGG,GGACGACAAGTGCAGG,GGATTGGCTCTGCAGG,GGCAAGGAGCTGCAGG,GGCACCTTACTGCAGG,GGCCATCCTCTGCAGG,GGCCGACAACTGCAGG,GGAGTCAAGCTGCAGG,GGCCTCTATCTGCAGG,GGCGACACACTGCAGG,GGCGGATTGCTGCAGG,GGCTAAGGTCTGCAGG,GGGAACAGGCTGCAGG,GGGACAGTGCTGCAGG,GGGAGTTAGCTGCAGG,GGGATGAATCTGCAGG,GGGCCAAGACTGCAGG ${i}_${c1}" >> 04_split_wells.${i}.${c1}.sh
echo "done creating script for" ${c1}
echo "Submitting job for 04_split_wells.${i}.${c1}.sh"
sbatch --mem MaxMemPerNode 04_split_wells.${i}.${c1}.sh
x=$(( $x + 1 ))
done

cd  $data_dir
done
