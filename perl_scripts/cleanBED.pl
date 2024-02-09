#!/usr/bin/perl
use strict;
use warnings;

my %chr;

open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t",$_);
	$chr{$col[0]} = $col[1];
}
close F;

my @file;

my $factor = $ARGV[1]/1000000;
print STDERR " - adjusted factor = $factor\n";
open G, $ARGV[2] or die;
while(<G>){
	chomp;
	my @col = split("\t",$_);
	if($col[2] > $chr{$col[0]}){
		$col[2] = $chr{$col[0]} - 1;
		if(($col[2] - $col[1]) < 1){
			next;
		}
	}
	if($col[1] >= $chr{$col[0]}){
		next;
	}
	my $adj = $col[3] / $factor;
	print "$col[0]\t$col[1]\t$col[2]\t$adj\n";
}
close G;

