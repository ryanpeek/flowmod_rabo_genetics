#!/usr/bin/perl
# Script written in 2014 by Mike Miller (micmiller@ucdavis.edu)
# Updated May 2019 by Shannon Joslin (sejoslin@ucdavis.edu)

# This script will split reads based on RAD barcodes.
# Usage: perl /path/to/BarcodeSplitPairedEndRAD.2019.pl <fastaR1> <fastaR3> <list,of,comma,separated,barcodes> <outputPrefix>


if ($#ARGV == 3) {
	$file1 = $ARGV[0];
	$file2 = $ARGV[1];
	$barcode = $ARGV[2];
	$prefix = $ARGV[3];
} else {
	die;	
}


@commas = split(/\,/, $barcode); # create an array of barcodes
$barcode_length = length($commas[0]); # capture length of barcode

print "Barcode length is $barcode_length\n";

$x=0;
# note: $#array gives index of last element
while ($x <= $#commas) {
        # make files names from hashes where key of hash refers to name of file
	$hash_r1{$commas[$x]} = $prefix . "_RA_" . $commas[$x] . ".fastq";
	$hash_r2{$commas[$x]} = $prefix . "_RB_" . $commas[$x] . ".fastq";
	# turn each element of the hash into a variable containing filenames
        $filename_r1 = $hash_r1{$commas[$x]};
        $filename_r2 = $hash_r2{$commas[$x]};
	# open one RA and RB file for each barcode
        open($filename_r1, ">$filename_r1") or die;
        open($filename_r2, ">$filename_r2") or die;
        $x++;
}

open(FILE1, "<$file1") or die;
open(FILE2, "<$file2") or die;

# while loop will store and read each line
#	and will continue until each line has been read
while (<FILE1>) {
	# grab R1 lines
        $f1a = $_;
        $f1b = <FILE1>;
        $f1c = <FILE1>;
        $f1d = <FILE1>;

	# grab R3 lines
        $f2a = <FILE2>;
        $f2b = <FILE2>;
        $f2c = <FILE2>;
        $f2d = <FILE2>;

	# pull barcodes from the sequence line (second fastq line 2) of file 1 and file 2
        $bc1 = substr($f1b,0,$barcode_length);
        $bc2 = substr($f2b,0,$barcode_length);

	# save string without barcodes
	# FILE1
	$f1b_2 = substr($f1b, $barcode_length, length($f1b));
	$f1d_2 = substr($f1d, $barcode_length, length($f1d));
	# FILE2
	$f2b_2 = substr($f2b, $barcode_length, length($f2b));
	$f2d_2 = substr($f2d, $barcode_length, length($f2d));


	# if the r1 hash contains something and r2 doesn't, then take the last 134 bps and save to variable name
        if ($hash_r1{$bc1} ne "" && $hash_r1{$bc2} eq "")  {

		# store file name as a variable
                $out1 = $hash_r1{$bc1};
                $out2 = $hash_r2{$bc1};
	
		# print shortened reads to files
                print $out1 $f1a . $f1b_2 . $f1c . $f1d_2;
                print $out2 $f2a . $f2b_2 . $f2c . $f2d_2;

        } elsif ($hash_r1{$bc1} eq "" && $hash_r1{$bc2} ne "")  {

                $out1 = $hash_r1{$bc2};
                $out2 = $hash_r2{$bc2};

                print $out1 $f2a . $f2b_2 . $f2c . $f2d_2;
                print $out2 $f1a . $f1b_2 . $f1c . $f1d_2;

        } elsif ($hash_r1{$bc1} ne "" && $hash_r1{$bc2} ne "")  {

                print "Double Barcode!\t$bc1\t$bc2\n";

        } elsif ($hash_r1{$bc1} eq "" && $hash_r1{$bc2} eq "") {

		print "Barcode not found!\t$bc1\t$bc2\n";
	}

}
close FILE1; close FILE2;

