#!/usr/bin/perl
use strict;
use warnings;
open my $out, ">", "Genes_replaced_in_trees.txt";
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	my @tarray=split/\s+/;
	
	my $test= `grep "$tarray[0]" good_pairs_james_maize.txt`;
	if($test){
		print $out "$tarray[0] replaced by $tarray[2]\n";
		`perl -pi -e "s/$tarray[0]/$tarray[2]/g" good_pairs_james_maize.txt`;
	}
}