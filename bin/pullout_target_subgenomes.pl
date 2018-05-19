#!/usr/bin/perl -w
use strict;

open my $file, "<", $ARGV[0];

open my $out, ">", "Target_gene_subgenomes.txt";
open my $out3, ">", "Single_subgenome_target.txt";
open my $out5, ">", "Single_subgenome_fourtargets.txt";
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
while(<$file>){
		chomp;
		my @tarray = split /\s+/;
		if($tarray[1] < 5){
			next;
		}
		if($tarray[6] =~ /Panicum/ || $tarray[6] =~ /Setaria/){
			next;
		}
		if($tarray[4] =~ /,/){
			my @temp1 = split(/,/, $tarray[4]);
			
			if($tarray[6] =~ /,/){
				my @temp2 = split(/,/, $tarray[6]);

				for my $t1 (@temp1){
					for my $t2 (@temp2){
						
						print $out "$t1\t$t2\n";
						$t2 =~ /^(.*?)\./;
						unless(exists $used{$t2}){
							$unique{$1}++;
							$used{$t2}=1;
						}
					}
				}
			}
			else{
				for my $t1 (@temp1){
					print $out "$t1\t$tarray[6]\n";
					print $out3 "$t1\t$tarray[6]\n";
					if($tarray[6] =~ /Vossia|Urely|Zmays|Hem|Coel/){
						if($t1 =~ /Zmays/){
							print $out5 "$t1\t$tarray[6]\n";
						}
					}
					unless($t1 =~ /Zmays/){
						next;
					}
					my $chrome=$chromosomes{$t1};
					$tarray[6] =~ /^(.*?)\./;
					my $species_n=$1;
					print $out6 "$species_n\t$chrome\t$t1\n";

					}
					
					$tarray[6] =~ /^(.*?)\./;
						unless(exists $used{$tarray[6]}){
							$unique{$1}++;
							$used{$tarray[6]}=1;
							$single_used{$1}++;
						}
				
			}
		}
		else{
			if($tarray[6] =~ /,/){
				my @temp2 = split(/,/, $tarray[6]);

				for my $t2 (@temp2){
						print $out "$tarray[4]\t$t2\n";
						if($tarray[6] =~ /Vossia|Urely|Zmays|Hem|Coel/){
							if($tarray[4] =~ /Zmays/){

								print $out5 "$tarray[4]\t$tarray[6]\n";
							}
							unless($tarray[4] =~ /Zmays/){
								next;
							}
						}
						my $chrome=$chromosomes{$tarray[4]};
						$t2 =~ /^(.*?)\./;
						my $species_n=$1;
							print $out6 "$species_n\t$chrome\t$tarray[4]\n";

					
						
						unless(exists $used{$t2}){
							$unique{$1}++;
							$used{$t2}=1;
						}
				}
				
			}
			else{
					print $out "$tarray[4]\t$tarray[6]\n";
					print $out3 "$tarray[4]\t$tarray[6]\n";
					if($tarray[6] =~ /Vossia|Urely|Zmays|Hem|Coel/){
						if($tarray[4] =~ /Zmays/){

							print $out5 "$tarray[4]\t$tarray[6]\n";
						}
					}
						unless($tarray[4] =~ /Zmays/){
							next;
						}
						my $chrome=$chromosomes{$tarray[4]};
						$tarray[6] =~ /^(.*?)\./;
						my $species_n=$1;
						print $out6 "$species_n\t$chrome\t$tarray[4]\n";
						
					}

					$tarray[6] =~ /^(.*?)\./;
					unless(exists $used{$tarray[6]}){
							$unique{$1}++;
							$used{$tarray[6]}=1;
							$single_used{$1}++;
						}
				}
			
		}

open my $out4, ">", "Unique_single_subgenome_contributions.txt";
for my $sp (sort keys %single_used){
		print $out4 "$sp\t$single_used{$sp}\n";
}

open my $out2, ">", "Unique_subgenome_contributions.txt";
for my $sp (sort keys %unique){
		print $out2 "$sp\t$unique{$sp}\n";
}
