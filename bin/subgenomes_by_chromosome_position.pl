#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0]; #read in chromosome file

my %gene2subgenome;

while(<$file>){
	chomp;
	my @tarray = split /\s+/;
	
	$tarray[2]=~ /Zmays.(.*?)$/;
	my $zea = $1;
	my $genus;
	if($tarray[0] =~ /-/){
		$tarray[0] =~ /(.*?)-/;
		$genus = $1;
	}
	else{
		$genus = $tarray[0];
	}
	$gene2subgenome{$zea}{$genus}=1;
}

open my $out, ">", "Putative_chromosome_ID_BSV$ARGV[2].txt";
print $out "Chromosome\tGene_ID\tStart\tPutative_Sub\n";
open my $gff, "<", $ARGV[1]; #gff file
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
			for my $id (sort {$a cmp $b} keys %{$gene2subgenome{$genetemp}}){
				push (@subs, $id);
			}
			my $sub = join(",", @subs);
			my $ctest=0;
			for my $ts (@subs){
				if($ts !~ /Urel/ && $ts !~ /Voss/){
					$ctest++;
				}
			}
			#if($ctest > 0){
			#	$sub = "Other";
			#} 
			print $out "$tarray[0]\t$genetemp\t$tarray[3]\t$sub\n";
		}
		else{
			print "Can't find $genetemp\n";
		}
	}
}
