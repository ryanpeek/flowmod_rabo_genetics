#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

while (<FILE>) {
	$line = $_;
	@tabs = split(/\t/,$line);
	$x = 2;
	$sum = 0;
	$num = 0;
	while ($x <= $#tabs) {
		$sum = $sum + $tabs[$x];
		$num++;
		$x++;
	}
	$mean = $sum/$num;
	print "$tabs[0]\t$tabs[1]\t$mean\n";
}

close(FILE);
