#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0];

open my $out6, ">", "Chromosome_by_subgenome.txt";

my %unique;
my %used;
my %chromosomes;
my %single_used;

open my $chromfile, "<", $ARGV[1];
while(<$chromfile>){
	chomp;
	my @tarray = split /\s+/;
	$chromosomes{"Zmays.".$tarray[1]}=$tarray[0];
}
print $out6 "Chromosome_P1\tParalog1_Name\tP1Sister_Taxa\tChromosome_P2\tParalog2_Name\tP2Sister_Taxa\tOutgroup_Sisters\n";
while(<$file>){
		chomp;
		my @tarray = split /\s+/;
		if($tarray[2] < 5){
			next;
		}
		
		my @p1_array = split(/,/, $tarray[6]);
		my @p2_array = split(/,/, $tarray[9]);
		my @out_taxaA = split(/,/, $tarray[10]);
		my %p1_taxa;
		my %p2_taxa;
		my %out_taxaH;
		
		my $temp1;
		my $temp2;
		my $temp0;
		for my $tax0 (@out_taxaA){
			
			if($tax0 =~ /-/){
				$tax0 =~ /^(.*?)-/;
				$temp0=$1;
			}
			else{
				$tax0 =~ /^(.*?)\./;
				$temp0=$1;
			}
			$out_taxaH{$temp0}=1;
		}
		for my $tax1 (@p1_array){
			if($tax1 eq $tarray[4]){
				
				next;
			}
			if($tax1 =~ /-/){
				$tax1 =~ /^(.*?)-/;
				$temp1=$1;
			}
			else{
				$tax1 =~ /^(.*?)\./;
				$temp1=$1;
			}
			
			
				$p1_taxa{$temp1}=1;
			
		}
		if(scalar @p1_array ==1 ){
				$p1_taxa{"SINGLE"}=1;
			}
		for my $tax2 (@p2_array){
			if($tax2 eq $tarray[7]){
				next;
			}
			if($tax2 =~ /-/){
				$tax2 =~ /^(.*?)-/;
				$temp2=$1;
			}
			else{
				$tax2 =~ /^(.*?)\./;
				$temp2=$1;
			}
			
			
				$p2_taxa{$temp2}=1;
			
		}
		if(scalar @p2_array ==1 ){
				$p2_taxa{"SINGLE"}=1;
			}
		my @outT;
		for my $out1k (sort {$a cmp $b} keys %out_taxaH){
			push(@outT, $out1k);
		}
		my $OUTtax = join (",", @outT);
		my @p1t;
		for my $tax1k (sort {$a cmp $b} keys %p1_taxa){
			push(@p1t, $tax1k);
		}
		my $ptax1 = join (",", @p1t);
		
		my @p2t;
		for my $tax2k (sort {$a cmp $b} keys %p2_taxa){
			push(@p2t, $tax2k);
		}
		my $ptax2 = join (",", @p2t);
		
		print $out6 "$chromosomes{$tarray[4]}\t$tarray[4]\t$ptax1\t$chromosomes{$tarray[7]}\t$tarray[7]\t$ptax2\t$OUTtax\n";
}
		