#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0];

open my $out, ">", "Target_gene_subgenomes_BSV$ARGV[3].txt";
open my $out3, ">", "Single_subgenome_target_BSV$ARGV[3].txt";
open my $out5, ">", "Single_subgenome_fourtargets_BSV$ARGV[3].txt";
open my $out6, ">", "Chromosome_by_subgenome_BSV$ARGV[3].txt";
open my $out7, ">", "Chromosome_by_subgenome_ParalogHIT_BSV$ARGV[3].txt";

my %unique;
my %used;
my %chromosomes;
my %single_used;

open my $para, "<", $ARGV[2];
my %paralogs;
while(<$para>){
	chomp;
	my @tarray = split /\s+/;
	$paralogs{$tarray[1]}=$tarray[0];
	$paralogs{$tarray[0]}=$tarray[1];
}

open my $chromfile, "<", $ARGV[1];
while(<$chromfile>){
	chomp;
	my @tarray = split /\s+/;
	$chromosomes{"Zmays.".$tarray[1]}=$tarray[0];
}
while(<$file>){
		chomp;
		my @tarray = split /\s+/;
		if($tarray[1] < $ARGV[3]){
			next;
		}
		#if($tarray[7] =~ /Panicum/ || $tarray[7] =~ /Setaria/){
		#	next;
		#}
		if($tarray[5] =~ /,/){
			my @temp1 = split(/,/, $tarray[5]);
			
			if($tarray[7] =~ /,/){
				my @temp2 = split(/,/, $tarray[7]);

				for my $t1 (@temp1){
					if(exists $paralogs{$t1}){
						for my $t2 (@temp2){
							my $chrome=$chromosomes{$t1};
							$t2 =~ /^(.*?)\./;
							my $species_n=$1;
							if($tarray[2] eq "Yes"){
								
								print $out7 "Zmays\t$chrome\t$t1\n";
							}
							else{
								print $out6 "$species_n\t$chrome\t$t1\n";
							}

						}
					}	
						
						
				}
				
			}
			else{
				for my $t1 (@temp1){
					print $out "$t1\t$tarray[6]\n";
					print $out3 "$t1\t$tarray[6]\n";
					if($tarray[7] =~ /Vossia|Urely|Zmays|Hem|Coel/){
						if($t1 =~ /Zmays/){
							print $out5 "$t1\t$tarray[6]\n";
						}
					}
					unless(exists $paralogs{$t1}){
						next;
					}
					my $chrome=$chromosomes{$t1};
					$tarray[7] =~ /^(.*?)\./;
					my $species_n=$1;
					if($tarray[2] eq "Yes"){
								print $out7 "Zmays\t$chrome\t$t1\n";
							}
							else{
								print $out6 "$species_n\t$chrome\t$t1\n";
							}
					

					}
					
					$tarray[7] =~ /^(.*?)\./;
						unless(exists $used{$tarray[7]}){
							$unique{$1}++;
							$used{$tarray[7]}=1;
							$single_used{$1}++;
						}
				
			}
		}
		else{
			unless(exists $paralogs{$tarray[5]}){
				next;
			}
			if($tarray[7] =~ /,/){
				my @temp2 = split(/,/, $tarray[7]);

				for my $t2 (@temp2){
						print $out "$tarray[5]\t$t2\n";
						if($tarray[7] =~ /Vossia|Urely|Zmays|Hem|Coel/){
							if($tarray[5] =~ /Zmays/){

								print $out5 "$tarray[5]\t$tarray[7]\n";
							}
							unless($tarray[5] =~ /Zmays/){
								next;
							}
						}
						my $chrome=$chromosomes{$tarray[5]};
						$t2 =~ /^(.*?)\./;
						my $species_n=$1;
						if($tarray[2] eq "Yes"){
								print $out7 "Zmays\t$chrome\t$tarray[5]\n";
							}
							else{
								print $out6 "$species_n\t$chrome\t$tarray[5]\n";
							}
							

					
						
						unless(exists $used{$t2}){
							$unique{$1}++;
							$used{$t2}=1;
						}
				}
				
			}
			else{
					print $out "$tarray[5]\t$tarray[7]\n";
					print $out3 "$tarray[5]\t$tarray[7]\n";
					if($tarray[6] =~ /Vossia|Urely|Zmays|Hem|Coel/){
						if($tarray[4] =~ /Zmays/){

							print $out5 "$tarray[5]\t$tarray[7]\n";
						}
					}
						unless($tarray[5] =~ /Zmays/){
							next;
						}
						my $chrome=$chromosomes{$tarray[5]};
						$tarray[7] =~ /^(.*?)\./;
						my $species_n=$1;
						if($tarray[2] eq "Yes"){
								print $out7 "Zmays\t$chrome\t$tarray[5]\n";
							}
							else{
								print $out6 "$species_n\t$chrome\t$tarray[5]\n";
							}
						
						
					}

					$tarray[7] =~ /^(.*?)\./;
					unless(exists $used{$tarray[7]}){
							$unique{$1}++;
							$used{$tarray[7]}=1;
							$single_used{$1}++;
						}
				}
			
		}

open my $out4, ">", "Unique_single_subgenome_contributions_BSV$ARGV[3].txt";
for my $sp (sort keys %single_used){
		print $out4 "$sp\t$single_used{$sp}\n";
}

open my $out2, ">", "Unique_subgenome_contributions_BSV$ARGV[3].txt";
for my $sp (sort keys %unique){
		print $out2 "$sp\t$unique{$sp}\n";
}
