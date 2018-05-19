#!/usr/bin/env perl -w
use strict;
use Bio::TreeIO;
open my $para, "<", $ARGV[3];
my $count_par=0;
my %used;
my %paralogs;

open my $out1, ">", "Orthogroups_with_paralogs.txt";
open my $out2, ">", "Paralogs_missing_mate.txt";

while(<$para>){
	chomp;
	my @tarray = split /\s+/;
	$paralogs{$tarray[1]}=$tarray[0];
	$paralogs{$tarray[0]}=$tarray[1];
}

my $dir = $ARGV[0];
my $sample;
if($ARGV[1] =~ /,/){
	my @temp = split(",", $ARGV[1]);
	$sample = join("|", @temp);
}
else{
	$sample = $ARGV[1];
}
opendir my $dh, $dir or die "Couldn't open '$dir' for reading: $!\n";
my @files = grep { !/^\.{1,2}$/ } readdir $dh;
closedir($dh);

my %tree_data;
my $count_tree=0;

for my $file (@files){
	$file = $dir . "/" . $file;

	my $treeio = new Bio::TreeIO(-format => "newick", -file => "$file", -internal_node_id => 'bootstrap');
	my %used_taxa;
	
	my $tree_switch=0;
	while( my $tree = $treeio->next_tree ) {
		if($tree_switch == 1){
			$count_tree++;
			$tree_switch=0;
		}
		
		my %genes;
		my @taxa = $tree->get_leaf_nodes;
		for my $tax (@taxa){
			my $taxon_id = $tax->id;
			if ($taxon_id !~ /Zmays/ ){
					next;
			}
			$genes{$taxon_id}=1;
		}
		
		for my $gtax (keys %genes){
			if(exists $paralogs{$gtax} && exists $genes{$paralogs{$gtax}}){
				if(exists $used_taxa{$gtax} || exists $used_taxa{$paralogs{$gtax}}){
					next;
				}
				else{
					print $out1 "$file\t$gtax\t$paralogs{$gtax}\n";
					$used_taxa{$gtax}=1;
					$used_taxa{$paralogs{$gtax}}=1;
				}
			}
			else{
				if(!exists $paralogs{$gtax}){
					print "$gtax is not a paralog\n";
				}
				else{
					print $out2 "$paralogs{$gtax} is the paralog of $gtax but not found in $file\n";
				}
			}
		}
	}
}
		
		