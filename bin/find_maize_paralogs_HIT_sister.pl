#!/usr/bin/env perl
use strict;
use warnings;

open my $file, "<", $ARGV[0];
open my $out, ">", "Target_Zmays_paralogs_HIT.txt";

my %para;
while (<$file>){
	chomp;
	$para{"Zmays.".$_}=1;
}

open my $file2, "<", $ARGV[1];
while(<$file2>){
	chomp;
	my @tarray = split/\s+/;
	
	if (exists $para{$tarray[0]} || exists $para{$tarray[1]}){
		print $out "$_\n";
	}
}