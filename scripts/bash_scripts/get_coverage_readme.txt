# To Calc Coverage/Depth

1. First index alignment/genome reference:
	- srun -t 120 bwa index final_contigs_300.fasta
2. Then make a file of baits and positions with samtools:
	- srun -t 1:00:00 --mem=5G bwa mem final_contigs_300.fasta baits_120_nanorana_only_8578.fa > baits.sam
3. Edit the file to trim start of lines: `grep -v "^@SQ" baits.sam > baits.makedepth`
4. Now make a new file that lists the position of the 60th read at each bait (could use any position, but midway good):
	- srun -t 6:00:00 MakeDepthBed_60.pl baits.makedepth > baits.depth60
5. Use this file to make a file where each column is a different individual at position 60 in this case:
	- srun -t 360 samtools depth -b baits.depth60 -f SOMM446/bamlists/rapture_2019_all_plates.bamlist > SOMM446/results_depth/baits.depth60_rapture_2019_all
6. Use this file to calculate the mean coverage of each bait at position 60 or mean coverage of each individual:
	- srun -t 360 code/07_CalcMeanCoverage.pl SOMM446/results_depth/baits.depth60_rapture_2019_all > SOMM446/results_depth/rapture_2019_all_site.meancoverage
	- srun -t 360 code/07_CalcMeanSampleCoverage.pl SOMM446/results_depth/baits.depth60_rapture_2019_all > SOMM446/results_depth/rapture_2019_all_sample.meancoverage
