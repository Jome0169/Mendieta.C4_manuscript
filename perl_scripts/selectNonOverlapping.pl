#!/usr/bin/perl
use strict;
use warnings;

my $chrom;

open F, $ARGV[0] or die;
while(<F>){
	chomp;
	my @col = split("\t",$_);
	if($_ =~ /B73/){
		$chrom = $col[0];
		$_ =~ s/B73V4_ctg/B73V4ctg/g;
		@col = split("\t",$_);
	}else{
		$chrom = $col[0];
	}
	if($col[3] !~ /\,/){
		my @score = split("_",$col[3]);
		print "$chrom\t$col[1]\t$col[2]\t$score[$#score]\n";
	}else{
		my @sites = split(",", $col[3]);
		my %hash;
		foreach(@sites){
			my @coord = split("_",$_);
			$hash{$_} = $coord[3];
		}
		my @sorted = sort {$hash{$b} <=> $hash{$a}} keys %hash;
		my @ncor = split("_",$sorted[0]);
		print "$chrom\t$ncor[1]\t$ncor[2]\t$ncor[3]\n";
		my @taken;
		push(@taken, $sorted[0]);
		for (my $i = 1; $i < @sorted; $i++){
			my @pos = split("_",$sorted[$i]);
			my $overlap = 0;
			for (my $j = 0; $j < @taken; $j++){
				my @comp = split("_",$taken[$j]);
				if($pos[1] >= $comp[1] && $pos[1] <= $comp[2]){
					$overlap++;
				}elsif($pos[2] >= $comp[1] && $pos[2] <= $comp[2]){
					$overlap++;
				}elsif($comp[1] <= $pos[1] && $comp[2] >= $pos[2]){
					$overlap++;
				}
			}
			if($overlap == 0){
				print "$chrom\t$pos[1]\t$pos[2]\t$pos[3]\n";
				push(@taken, $sorted[$i]);
			}
		}
	}
}
close F;
