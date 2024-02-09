#!/usr/bin/perl
use strict;
use warnings;

my $total = 0;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t",$_);
	$total = $total + $col[4];
}
close F;

my $scalefactor = $total / 1000000;
open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t",$_);
	my $score = $col[4] / $scalefactor;
	print "$col[0]\t$col[1]\t$col[2]\t$score\n";
}
close F;

