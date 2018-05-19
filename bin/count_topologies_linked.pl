#!/usr/bin/perl
use strict;
use warnings;


my %counts;
open my $file, "<", $ARGV[0];

while(<$file>){
	chomp;
	next if /Chromosome/;
	my @tarray=split/\s+/;
	
	if(exists $counts{$tarray[3]} || exists $counts{$tarray[8]}){
		
		if(exists $counts{$tarray[3]}{$tarray[8]}){
			$counts{$tarray[3]}{$tarray[8]}++;
		}
		elsif(!exists $counts{$tarray[3]}{$tarray[8]}){
			$counts{$tarray[3]}{$tarray[8]}++;
		}
		elsif(exists $counts{$tarray[8]}{$tarray[3]}){
			$counts{$tarray[8]}{$tarray[3]}++;
		}
		else{
			$counts{$tarray[8]}{$tarray[3]}++;
		}
	}
	else{
		$counts{$tarray[3]}{$tarray[8]}++;
	}
}

open my $out, ">", "Results_for_BSV$ARGV[1].txt";
for my $set (keys %counts){
	for my $set2 (keys %{$counts{$set}}){
		print $out "$set\t$set2\t$counts{$set}{$set2}\n";
	}
}
		
	
	