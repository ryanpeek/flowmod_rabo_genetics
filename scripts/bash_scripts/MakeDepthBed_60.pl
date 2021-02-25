#!/usr/bin/perl

$file = $ARGV[0];

open(FILE, "<$file") or die;

	while (<FILE>) {
		@tabs = split(/\t/,$_);
		if ($tabs[1] == 16) {
			$start = $tabs[3]-14;
			$stop = $start + 133;
			$start2 = $stop-60;
			$stop2 = $start2+1;
		} elsif ($tabs[1] == 0) {
			$start = $tabs[3];
			$stop = $start + 133;
			$start2 = $start+59;
			$stop2 = $start2+1;
		}
		$chr = $tabs[2];
		if ($start2 < 1) { $start2=1; }
		if ($stop2 < 1) { $stop2=1; }
		print "$chr\t$start2\t$stop2\n";
	}

close(FILE);
