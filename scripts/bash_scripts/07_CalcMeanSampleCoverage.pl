#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

$l=0;
while (<FILE>) {
	$line = $_; chomp($line);
	@tabs = split(/\t/,$line);
	$x = 2;
	while ($x <= $#tabs) {
		$sum[$x-2] = $sum[$x-2] + $tabs[$x];
		$x++
	}
	$l++;
}
close(FILE);

$x=0;
while ($x <= $#sum) {
	$mean = $sum[$x]/$l;
	print "$x\t$mean\n";
	$x++;
}

