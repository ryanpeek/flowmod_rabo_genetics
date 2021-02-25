#!/bin/bash -l

#SBATCH --mail-user=rapeek@ucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH -J splitSEQ
#SBATCH -e slurms/05b_subsample.%j.err
#SBATCH -o slurms/05b_subsample.%j.out
#SBATCH -c 20
#SBATCH -p high
#SBATCH --time=2-20:00:00

set -e # exits upon failing command
set -v # verbose -- all lines
set -x # trace of all commands after expansion before execution#SBATCH -J subsample

# run with "sbatch 05b_run_subsample.sh bamlist 30k


list=$1
num=$2

wc=$(wc -l ${list} | awk '{print $1}')

x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1, $2}')   
	set -- $var
	c1=$1
	c2=$2

	count=$(samtools view -c ${c1})

        if [ $num -le $count ]
        then
                frac=$(bc -l <<< $num/$count)
                samtools view -bs $frac ${c1} > ${c1}_${num}.bam
        fi

	x=$(( $x + 1 ))

done


