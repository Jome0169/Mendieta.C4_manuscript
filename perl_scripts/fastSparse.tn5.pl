#!/usr/bin/perl
use strict;
use warnings;

my %hash;
my $its = 0;

open F, $ARGV[0] or die;
while(<F>){
	$its++;
	if(($its % 1000000)==0){
		print STDERR " - read $its reads ... \n";
	}
	chomp;
	my @col = split("\t",$_);
	$col[0] = 'chr' . $col[0];
	if($col[0] =~ /B73V/){
		$col[0] =~ s/chrB73V4_ctg/chrB73V4ctg/g;
	}
	my $id = join("_", $col[0],$col[6],$col[7]);
	if(exists $hash{$col[3]}{$id}){
		next;
	}else{
		print "$id\t$col[3]\t1\n";
		$hash{$col[3]}{$id} = 1;
	}
}
close F;
