#!/usr/bin/env perl -w
use strict;
use Bio::TreeIO;
open my $para, "<", $ARGV[3];
my $count_par=0;
my %used;
my %paralogs;
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
my %missing_pars;
my %split_pars;
for my $file (@files){
	$file = $dir . "/" . $file;

	my $treeio = new Bio::TreeIO(-format => "newick", -file => "$file", -internal_node_id => 'bootstrap');
	my %used_taxa;
	
	$missing_pars{$file}=0;
	while( my $tree = $treeio->next_tree ) {
		
		
			
		my @taxa = $tree->get_leaf_nodes;
		my %temp_taxa;
		for my $taxon (@taxa){
			$temp_taxa{$taxon->id}=1;
		}
			
		for my $taxon (@taxa){
			my $ancestor;
			my $bs;
			my $taxon_id = $taxon->id;
			
			if ($taxon_id !~ /Zmays/ ){
					next;
			}
			unless(exists $paralogs{$taxon_id}){
				next;
			}
			if(exists $used{$taxon_id}){
				next;
			}
			$missing_pars{$file}=1;
			if(exists $temp_taxa{$paralogs{$taxon_id}}){
				$split_pars{$file}+=0;
				
			}
			else{
				$split_pars{$file}+=1;
				$missing_pars{$file}+=0;
			}
			
			$count_par++;
			my $descendents;

			my (@descendents,@old_descendents);
			
			$ancestor = $taxon->ancestor;
			$bs = $ancestor->bootstrap;
			if(!$bs){
				$bs = 0;
			}
			my @temp_descendents = $ancestor->get_all_Descendents;
			
			for my $temp_descend (@temp_descendents){
				if($temp_descend->is_Leaf){
					$temp_descend=$temp_descend->id;
					push(@descendents, $temp_descend);
				}
			}	

			
		
			my $count = grep (/$sample/, @descendents);
			#print "$taxon_id\n";
			my $count_p = grep (/$paralogs{$taxon_id}/, @descendents);
			
			#Until our set of tips includes taxa outside of the focal taxon, we keep moving one node back in the tree.
			while ($count == @descendents && $count_p == 0){
				
				($descendents, $ancestor, $bs) = get_descendents(\@descendents, $ancestor);
				@descendents=@$descendents;
				$count = grep (/$sample/, @descendents);
				$count_p = grep (/$paralogs{$taxon_id}/, @descendents);
			}
			my $paralog_hit="No";
			if($count_p > 0){
				$paralog_hit="Yes";
			}
			#Getting the two child clades from our current ancestor.  Looking for components to identify if the focal taxon genes are more closely related to the alternative subgenome.
			my @subclade_ancestors;
			for my $child ($ancestor->each_Descendent){
				push(@subclade_ancestors, $child);
			}
			
			#Identifying tips for each clade.
			my @tips1;
			my @tips2;
			if($subclade_ancestors[0]->is_Leaf){
				my $single=$subclade_ancestors[0]->id;
				push(@tips1, $single);
				
			}
			else{
				my @clade1 = $subclade_ancestors[0]->get_all_Descendents;
			
				for my $temp_descend (@clade1){
					if($temp_descend->is_Leaf){
						$temp_descend=$temp_descend->id;
						push(@tips1, $temp_descend);
					}
				}
			}	
			if($subclade_ancestors[1]->is_Leaf){
				my $single=$subclade_ancestors[1]->id;
				push(@tips2, $single);

			}
			else{
				my @clade2 = $subclade_ancestors[1]->get_all_Descendents;
			
				for my $temp_descend (@clade2){
					if($temp_descend->is_Leaf){
						$temp_descend=$temp_descend->id;
						push(@tips2, $temp_descend);
					}
				}
			}	
			my $bs1 = $subclade_ancestors[0]->bootstrap;
			my $bs2 = $subclade_ancestors[1]->bootstrap;
			if (@tips1 == 1){
				$bs1 = "NA";
			} 
			if(@tips2 == 1) {
				$bs2 = "NA";
			}
			#Identifying if clades have target taxon.
			my $count1 = grep (/Zmays/, @tips1);
			my $count2 = grep (/Zmays/, @tips2);

			my $target_found="No";
			if($count1 > 0 && $count2 > 0){
				$target_found = "No";
			}
			
			#Determine which clade has the original tip in it.
			my $find_orig1 = grep (/$taxon_id/, @tips1);
			my $find_orig2 = grep (/$taxon_id/, @tips2);
			$used{$taxon_id}=1;

			if ($find_orig1 > 0){
				if(grep (/$paralogs{$taxon_id}/, @tips2) > 0){
					#$used{$paralogs{$taxon_id}}=1;
					$target_found = "Yes";
				}
				for my $tip (@tips1){
					push (@{$tree_data{$file}{$ancestor}{"Target"}{"Species"}}, $tip);
				}
				$tree_data{$file}{$ancestor}{"Target"}{"BS"}=$bs1;
				$tree_data{$file}{$ancestor}{"BS"}=$bs;
				$tree_data{$file}{$ancestor}{"BOTH"}=$target_found;
				$tree_data{$file}{$ancestor}{"Paralog"}=$paralog_hit;

				for my $tip (@tips2){
					push (@{$tree_data{$file}{$ancestor}{"Other"}{"Species"}}, $tip);
				}
				$tree_data{$file}{$ancestor}{"Other"}{"BS"}=$bs2;
				

			}
			else{
				for my $tip (@tips2){
					push (@{$tree_data{$file}{$ancestor}{"Target"}{"Species"}}, $tip);
				}
				if(grep (/$paralogs{$taxon_id}/, @tips1) > 0){
					$target_found = "Yes";
					#$used{$paralogs{$taxon_id}}=1;
				}
				$tree_data{$file}{$ancestor}{"Target"}{"BS"}=$bs2;
				$tree_data{$file}{$ancestor}{"BS"}=$bs;
				$tree_data{$file}{$ancestor}{"BOTH"}=$target_found;
				$tree_data{$file}{$ancestor}{"Paralog"}=$paralog_hit;
				

				for my $tip (@tips1){
					push (@{$tree_data{$file}{$ancestor}{"Other"}{"Species"}}, $tip);
				}
				$tree_data{$file}{$ancestor}{"Other"}{"BS"}=$bs1;

			}
		}
}
}

open my $out, ">", "Subgenome_results_for_" . $ARGV[2] . ".txt";
print $out "Tree\tNode_BS\tTarget_Both_Subclades\tTarget_Subclade_BS\tTarget_Subclade_Species\tTarget_Subclade_BS\tParalog_Hit_First\tOther_Subclade_BS\tOther_Subclade_Species\n";
for my $tree_file (sort keys %tree_data){
	for my $anc (keys %{$tree_data{$tree_file}}){
	my $sc1 = join (",", @{$tree_data{$tree_file}{$anc}{Target}{Species}});
	my $sc2 = join (",", @{$tree_data{$tree_file}{$anc}{Other}{Species}});

	print $out "$tree_file\t$tree_data{$tree_file}{$anc}{BS}\t$tree_data{$tree_file}{$anc}{BOTH}\t$tree_data{$tree_file}{$anc}{Target}{BS}\t$tree_data{$tree_file}{$anc}{Paralog}\t$sc1\t$tree_data{$tree_file}{$anc}{Other}{BS}\t$sc2\n";
	}
}
for my $tfname (keys %missing_pars){
	$count_tree+=$missing_pars{$tfname};
}

my $split_partree=0;
for my $sfname (keys %split_pars){
	$split_partree+=$split_pars{$sfname};
}

print "Total paralogs found: $count_par\nTrees with paralogs: $count_tree\nTrees with Split Paralogs:$split_partree\n";

sub get_descendents {
	#takes array of descendents and ancestor
	my $array_ref = $_[0];
	my @old_descendents= $array_ref;
	my $ancestor = $_[1]->ancestor;
	my $bs = $ancestor->bootstrap;
	if(!$bs){
			$bs = 0;
	}
	my @temp_descendents = $ancestor->get_all_Descendents;
	my @descendents=();
	for my $temp_descend (@temp_descendents){
		if($temp_descend->is_Leaf){
			$temp_descend=$temp_descend->id;
			push(@descendents, $temp_descend);
		}
	}	
	return (\@descendents, $ancestor, $bs);

}
