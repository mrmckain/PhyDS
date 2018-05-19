#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0]; #read in chromosome file

my %gene2subgenome;

while(<$file>){
	chomp;
	my @tarray = split /\s+/;
	
	$tarray[2]=~ /Zmays.(.*?)$/;
	my $zea = $1;
	$gene2subgenome{$zea}{$tarray[0]}=1;
}
my %paralogs;
open my $paralogfile, "<", $ARGV[1];
while(<$paralogfile>){
	chomp;
	my @tarray= split /\s+/;
	$paralogs{$tarray[0]}=$tarray[1];
	$paralogs{$tarray[1]}=$tarray[0];
}
open my $out, ">", "Putative_chromosome_ID.txt";
print $out "Chromosome\tGene_ID\tStart\tPutative_Sub\n";
open my $gff, "<", $ARGV[1]; #gff file
my %used;
while(<$gff>){
	chomp;
	if(/\#\#/){
		next;
	}
	
	my @tarray = split /\s+/;
	if($tarray[2] eq "gene"){
		$tarray[8] =~ /Name\=(.*?)$/;
		my $genetemp = $1;
		if(exists $gene2subgenome{$genetemp}){
			my @subs;
			for my $id (keys %{$gene2subgenome{$genetemp}}){
				push (@subs, $id);
			}
			my $sub = join(",", @subs);
			my $ctest=0;
			for my $ts (@subs){
				if($ts !~ /Urel/ && $ts !~ /Voss/){
					$ctest++;
				}
			}
			if($ctest > 0){
				$sub = "Other";
			} 
			print $out "$tarray[0]\t$genetemp\t$tarray[3]\t$sub";
			if(exists $gene2subenome{$paralogs{$genetemp}}){
				my @tsubs;
                        	for my $id (keys %{$gene2subgenome{$paralogs{$genetemp}}}){
                                	push (@tsubs, $id);
                        	}
                        	my $tsub = join(",", @subs);
                        	my $ctest=0;
                        	for my $ts (@tsubs){
                                	if($ts !~ /Urel/ && $ts !~ /Voss/){
                                        	$ctest++;
                                	}
                        	}
                        	if($ctest > 0){
                                	$sub = "Other";
                        	}
				print $out "
			}
		else{
			print "Can't find $genetemp\n";
		}
	}
}
