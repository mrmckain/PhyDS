#!/usr/bin/env perl
use warnings;
use strict;
use Bio::TreeIO;
use Getopt::Long;

my $paralog_file;
my $tree_dir;
my $help;
my $version;
my $tree_type="best";
my $ignore_taxa;
my $name = "PhyDS";
my $bootstrap_value = 0;
GetOptions('help|?' => \$help,'version' => \$version, "paralog_file=s" => \$paralog_file, "trees=s" => \$tree_dir, "bootstrap=i" => \$bootstrap_value, "ignore=s" => \$ignore_taxa, "name=s" => \$name, "tree_type=s" => \$tree_type)  or pod2usage( { -message => "ERROR: Invalid parameter." } );
my %final_counts;

open my $para, "<", $paralog_file;
my $target_id;
my $count_par=0;
my %used;
my %paralogs;
my %sister_taxa;
while(<$para>){
	chomp;
	my @tarray = split /\s+/;

	$paralogs{$tarray[1]}=$tarray[0];
	$paralogs{$tarray[0]}=$tarray[1];
	$tarray[0] =~ /^(.*?)-.+/;
	$target_id = $1;
}

my $dir = $tree_dir;
my $sample;
if($ignore_taxa =~ /,/){
	my @temp = split(",", $ignore_taxa);
	$sample = join("|", @temp);
}
else{
	$sample = $ignore_taxa;
}
opendir my $dh, $dir or die "Couldn't open '$dir' for reading: $!\n";
my @files = grep { !/^\.{1,2}$/ } readdir $dh;
closedir($dh);


if($tree_type eq "best"){
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
		
		
			
			my @taxa = $tree->get_leaf_nodes; #Get leaves
			my %temp_taxa;
			for my $taxon (@taxa){
				$temp_taxa{$taxon->id}=1; #get human readable ids of leaves
			}
			
			for my $taxon (@taxa){
				my $ancestor;
				
				my $bs;
				my $taxon_id = $taxon->id;
			
				#check to see if the leaf in question is a paralog or been used. Should check to make sure this isn't a problem with paralogs found in multiple events.
				unless(exists $paralogs{$taxon_id}){
					next;
				}
				if(exists $used{$taxon_id}){
					next;
				}
				$missing_pars{$file}=1;
				#checking if paralog mate present
				if(exists $temp_taxa{$paralogs{$taxon_id}}){
					$split_pars{$file}+=0;
				
				}
				else{
					$split_pars{$file}+=1;
					print "$taxon_id\n";
					$missing_pars{$file}+=0;
					next;
				}
			
				$count_par++;
				my $descendents;
				#initializing the descendents array for the MRCA of the paralog in question
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
				#The idea here is that we want one of two things: 1) we find a taxon that does not share the WGD event (i.e. not from our ignore list) or 2) is our putative paralog being the subgenomes are more closely related to each other.
				while ($count == @descendents && $count_p == 0){
				
					($descendents, $ancestor, $bs) = get_descendents(\@descendents, $ancestor);
					@descendents=@$descendents;
					#check if the taxa we are looking at are from the ignore list and count them. We are trying to determine if we have found taxa we don't want to ignore.
					$count = grep (/$sample/, @descendents);
					#count occurrences of putative paralog in the descendents
					$count_p = grep (/$paralogs{$taxon_id}/, @descendents);
				}
				my $paralog_hit="No";
				#check if both members of the paralog pair in the descendents
				my @sister_clade;

				if($count_p > 0 && $count == @descendents){
					$paralog_hit="Yes";
					#Get the sister clade if we have a clade of just the paralogs and ignore taxa.
					my($a_descendents, $a_ancestor, $a_bs) = get_descendents(\@descendents, $ancestor);
					my @a_descendents=@$a_descendents;
					my %seen;
					@seen{@descendents} = 1;
					my @solo_sister_des;
                    for my $item (@a_descendents) {
                     push(@sister_clade, $item) unless exists $seen{$item};
                    }
                   

					@sister_clade = sort {$a cmp $b} @sister_clade;
					if(!$a_bs){
						$a_bs = 0;
					}
					$sister_taxa{$file}{$a_ancestor}{"BS"}=$a_bs;

					$sister_taxa{$file}{$a_ancestor}{"TAXA"}=\@sister_clade;

				}
				if($count_p > 0 && $count != @descendents){
					$paralog_hit="Yes_other";
}
				#Getting the two child clades from our current ancestor.  Looking for components to identify if the focal taxon genes are more closely related to the alternative subgenome.
				my @subclade_ancestors;
				for my $child ($ancestor->each_Descendent){
					push(@subclade_ancestors, $child);
				}
			
			#Identifying tips for each clade.
				my @tips1;
				my @tips2;
				
				#first checking if a descendent of the current ancestor is a tip. if true, we store that as a single entry for that subclade
				if($subclade_ancestors[0]->is_Leaf){
					my $single=$subclade_ancestors[0]->id;
					push(@tips1, $single);
				
				}
				else{
					my @clade1 = $subclade_ancestors[0]->get_all_Descendents;
					#if the descendent of the ancestor is an internal node, we get all of the descedents are proceed through to get the tips
					for my $temp_descend (@clade1){
						if($temp_descend->is_Leaf){
							$temp_descend=$temp_descend->id;
							push(@tips1, $temp_descend);
						}
					}
				}	

				#repeat above for the other descendent of the ancestor
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
				#get bootstraps for the two clades descendent from ancestor node
				my $bs1 = $subclade_ancestors[0]->bootstrap;
				my $bs2 = $subclade_ancestors[1]->bootstrap;
				if (@tips1 == 1){
					$bs1 = "NA";
				} 
				if(@tips2 == 1) {
					$bs2 = "NA";
				}
			#Identifying if clades have target taxon. The target is the paralog we are focused on.
			#Check this. This code might be outdated since the recent update. mrm-2/17/2019
				my $count1 = grep (/$target_id/, @tips1);
				my $count2 = grep (/$target_id/, @tips2);

				my $target_found="No";
				if($count1 > 0 && $count2 > 0){
					$target_found = "No";
				}
			
			#Determine which clade has the original tip in it.
				my $find_orig1 = grep (/$taxon_id/, @tips1);
				my $find_orig2 = grep (/$taxon_id/, @tips2);
				$used{$taxon_id}=1;

				if ($find_orig1 > 0){
					if(grep (/$paralogs{$taxon_id}/, @tips2) > 0 && $count == @descendents){
					$used{$paralogs{$taxon_id}}=1;
						$target_found = "Yes";
					}
					if(grep (/$paralogs{$taxon_id}/, @tips2) > 0 && $count != @descendents){
					#$used{$paralogs{$taxon_id}}=1;
						$target_found = "Yes_other";
					}
					for my $tip (@tips1){
						push (@{$tree_data{$file}{$ancestor}{"Target"}{"Species"}}, $tip);
					}
					$tree_data{$file}{$ancestor}{"Target"}{"BS"}=$bs1;
					$tree_data{$file}{$ancestor}{"BS"}=$bs;
					$tree_data{$file}{$ancestor}{"BOTH"}=$target_found;
					$tree_data{$file}{$ancestor}{"Paralog"}=$taxon_id;

					for my $tip (@tips2){
						push (@{$tree_data{$file}{$ancestor}{"Other"}{"Species"}}, $tip);
					}
					$tree_data{$file}{$ancestor}{"Other"}{"BS"}=$bs2;
				

				}
				else{
					for my $tip (@tips2){
						push (@{$tree_data{$file}{$ancestor}{"Target"}{"Species"}}, $tip);
					}
					if(grep (/$paralogs{$taxon_id}/, @tips1) > 0 && $count == @descendents){
						$target_found = "Yes";
					$used{$paralogs{$taxon_id}}=1;
					}
					if(grep (/$paralogs{$taxon_id}/, @tips1) > 0 && $count != @descendents){
					#$used{$paralogs{$taxon_id}}=1;
						$target_found = "Yes_other";
					}
					$tree_data{$file}{$ancestor}{"Target"}{"BS"}=$bs2;
					$tree_data{$file}{$ancestor}{"BS"}=$bs;
					$tree_data{$file}{$ancestor}{"BOTH"}=$target_found;
					$tree_data{$file}{$ancestor}{"Paralog"}=$taxon_id;
				

					for my $tip (@tips1){
						push (@{$tree_data{$file}{$ancestor}{"Other"}{"Species"}}, $tip);
					}
					$tree_data{$file}{$ancestor}{"Other"}{"BS"}=$bs1;

				}
			}
		}
	}



	open my $out, ">", $name .  "_Subgenome_Results.txt";
	print $out "Tree\tBSV_MRCA\tBoth_Putative_Paralogs_Identified_in_Separate_Subtrees\tTarget_Paralog\tBSV_Subtree_with_Target_Paralog\tTaxa_in_Subtree_with_Target_Paralog\tBSV_Non-target_Subtree\tTaxa_in_Subtree_without_Target_Paralog\n";
	
	my (@sc1a, @sc2a);
	for my $tree_file (sort keys %tree_data){
		for my $anc (keys %{$tree_data{$tree_file}}){
			
		my $sc1 = join (",", @{$tree_data{$tree_file}{$anc}{Target}{Species}});
		my $sc2 = join (",", @{$tree_data{$tree_file}{$anc}{Other}{Species}});
		if($tree_data{$tree_file}{$anc}{BS} >= $bootstrap_value){ 
			if($tree_data{$tree_file}{$anc}{BOTH} eq "Yes"){
		@sc1a = sort {$a cmp $b} @{$tree_data{$tree_file}{$anc}{Target}{Species}};
		@sc2a = sort {$a cmp $b} @{$tree_data{$tree_file}{$anc}{Other}{Species}};
			my @temp_sc1 = @sc1a;
					my @temp_sc2 = @sc2a;

					@temp_sc1 = sort {$a cmp $b} @temp_sc1;
					@temp_sc2 = sort {$a cmp $b} @temp_sc2;

					my %temp_ids1;
					my %temp_ids2;
					for my $sp (@temp_sc1){
							$sp =~ /^(.*?)-/;
							$temp_ids1{$1}=1;

					}
					for my $sp (@temp_sc2){
							$sp =~ /^(.*?)-/;
							$temp_ids2{$1}=1;

					}

					my %names_temp1;
					my %names_temp2;
					for my $k (keys %temp_ids1){
							if($k =~ /-/){
									$k =~ /^(.*?)-/;
									$k = $1;
							}
							$names_temp1{$k}=1;
					}
					for my $k (keys %temp_ids2){
							if($k =~ /-/){
									$k =~ /^(.*?)-/;
									$k = $1;
							}
							$names_temp2{$k}=1;
					}

					my @names_temp1 = sort (keys %names_temp1);
					my @names_temp2 = sort (keys %names_temp2);
					my $tsp1 = join (",", @names_temp1);
					$final_counts{$tsp1}++;

					my $tsp2 = join (",", @names_temp2);
					$final_counts{$tsp2}++;
				}
				elsif($tree_data{$tree_file}{$anc}{BOTH} eq "Yes_other"){
					@sc1a = sort {$a cmp $b} @{$tree_data{$tree_file}{$anc}{Target}{Species}};
						my @fused = sort {$a cmp $b} @sc1a;

						my %temp_ids1;
						for my $sp (@fused){
							$sp =~ /^(.*?)-/;
							$temp_ids1{$1}=1;

					}
					my %names_temp1;
					for my $k (keys %temp_ids1){
							if($k =~ /-/){
									$k =~ /^(.*?)-/;
									$k = $1;
							}
							$names_temp1{$k}=1;
					}
					my @names_temp1 = sort (keys %names_temp1);
					my $tsp1 = join (",", @names_temp1);
					$final_counts{$tsp1}++;
				}
				else{
						@sc1a = sort {$a cmp $b} @{$tree_data{$tree_file}{$anc}{Target}{Species}};
						@sc2a = sort {$a cmp $b} @{$tree_data{$tree_file}{$anc}{Other}{Species}};
						my @fused = sort {$a cmp $b} (@sc1a,@sc2a);

						my %temp_ids1;
						for my $sp (@fused){
							$sp =~ /^(.*?)-/;
							$temp_ids1{$1}=1;

					}
					my %names_temp1;
					for my $k (keys %temp_ids1){
							if($k =~ /-/){
									$k =~ /^(.*?)-/;
									$k = $1;
							}
							$names_temp1{$k}=1;
					}
					my @names_temp1 = sort (keys %names_temp1);
					my $tsp1 = join (",", @names_temp1);
					$final_counts{$tsp1}++;
				}
		print $out "$tree_file\t$tree_data{$tree_file}{$anc}{BS}\t$tree_data{$tree_file}{$anc}{BOTH}\t$tree_data{$tree_file}{$anc}{Paralog}\t$tree_data{$tree_file}{$anc}{Target}{BS}\t$sc1\t$tree_data{$tree_file}{$anc}{Other}{BS}\t$sc2\n";
		}}
	}
	for my $tfname (keys %missing_pars){
		$count_tree+=$missing_pars{$tfname};
	}

	my $split_partree=0;
	for my $sfname (keys %split_pars){
		$split_partree+=$split_pars{$sfname};
	}
	open my $out_f, ">", "Subgenome_Counts_BSV-Cutoff_" . $bootstrap_value . "_" . $name . ".txt";
	for my $kid (sort {$a cmp $b} keys %final_counts){
			print $out_f "$kid\t$final_counts{$kid}\n";
	}
	print "Total paralogs found: $count_par\nTrees with paralogs: $count_tree\nTrees with Split Paralogs:$split_partree\n";

	open my $out_sister, ">", $name. "Taxa_Sister_to_ParalogsOnly_Clade.txt";
	my %sis_count;
	for my $tree (sort {$a cmp $b} keys %sister_taxa){
		for my $anc (keys %{$sister_taxa{$tree}}){
			if($sister_taxa{$tree}{$anc}{"BS"} >= 1){
				my %names_temp1;
				for my $k (@{$sister_taxa{$tree}{$anc}{"TAXA"}}){
							if($k =~ /-/){
									$k =~ /^(.*?)-/;
									$k = $1;
							}
							$names_temp1{$k}=1;
					}
				my @names_temp1 = sort (keys %names_temp1);
					my $tsp1 = join (",", @names_temp1);
				$sis_count{$tsp1}++;
			}
		}
	}

	for my $scount (sort {$a cmp $b} keys %sis_count){
			print $out_sister "$scount\t$sis_count{$scount}\n";
	}
}

if($tree_type eq "bootstrap"){
	my %tree_data;
	my $count_tree=0;
	my %missing_pars;
	my %split_pars;
	for my $file (@files){
		$file = $dir . "/" . $file;

		my $treeio = new Bio::TreeIO(-format => "newick", -file => "$file");
		
		my %used_taxa;
	
		
		my $tree_counter=0;
		while( my $tree = $treeio->next_tree ) {
			$missing_pars{$file}=0;
			$tree_counter++;
			#print "$file\n";
			my @taxa = $tree->get_leaf_nodes;
			my %temp_taxa;
			for my $taxon (@taxa){
				$temp_taxa{$taxon->id}=1;
			}
			
			for my $taxon (@taxa){
				my $ancestor;
			
				my $taxon_id = $taxon->id;
			
				#if ($taxon_id !~ /Zmays/ ){
				#	next;
				#}
				unless(exists $paralogs{$taxon_id}){
					next;
				}
				#if(exists $used{$taxon_id}){
				#	next;
				#}
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
						my $bs;
					($descendents, $ancestor) = get_descendents_boot(\@descendents, $ancestor);
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
				#my $bs1 = $subclade_ancestors[0]->bootstrap;
				#my $bs2 = $subclade_ancestors[1]->bootstrap;
				#if (@tips1 == 1){
				#	$bs1 = "NA";
				#} 
				#if(@tips2 == 1) {
				#	$bs2 = "NA";
				#}
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
						push (@{$tree_data{$file}{$tree_counter}{$ancestor}{"Target"}{"Species"}}, $tip);
					}
					$tree_data{$file}{$tree_counter}{$ancestor}{"Target"}{"BS"}++;
					$tree_data{$file}{$tree_counter}{$ancestor}{"BS"}++;
					$tree_data{$file}{$tree_counter}{$ancestor}{"BOTH"}=$target_found;
					$tree_data{$file}{$tree_counter}{$ancestor}{"Paralog"}=$paralog_hit;

					for my $tip (@tips2){
						push (@{$tree_data{$file}{$tree_counter}{$ancestor}{"Other"}{"Species"}}, $tip);
					}
					$tree_data{$file}{$tree_counter}{$ancestor}{"Other"}{"BS"}++;
				

				}
				else{
					for my $tip (@tips2){
						push (@{$tree_data{$file}{$tree_counter}{$ancestor}{"Target"}{"Species"}}, $tip);
					}
					if(grep (/$paralogs{$taxon_id}/, @tips1) > 0){
						$target_found = "Yes";
					#$used{$paralogs{$taxon_id}}=1;
					}
					$tree_data{$file}{$tree_counter}{$ancestor}{"Target"}{"BS"}++;
					$tree_data{$file}{$tree_counter}{$ancestor}{"BS"}++;
					$tree_data{$file}{$tree_counter}{$ancestor}{"BOTH"}=$target_found;
					$tree_data{$file}{$tree_counter}{$ancestor}{"Paralog"}=$paralog_hit;
				

					for my $tip (@tips1){
						push (@{$tree_data{$file}{$tree_counter}{$ancestor}{"Other"}{"Species"}}, $tip);
					}
					$tree_data{$file}{$tree_counter}{$ancestor}{"Other"}{"BS"}++;

				}
			
			}
		}
	}


	open my $out, ">", $name . "_Subgenome_Results.txt";
		print $out "Tree\tNode_Count\tTarget_Both_Subclades\tTarget_Subclade_Count\tTarget_Subclade_Species\tTarget_Subclade_Count\tParalog_Hit_First\tOther_Subclade_BS\tOther_Subclade_Species\n";

	for my $tree_file (sort keys %tree_data){
		for my $t_count (sort {$a <=> $b} keys %{$tree_data{$tree_file}}){
				print "$t_count\n";
		for my $anc (keys %{$tree_data{$tree_file}{$t_count}}){
			my @sc1a = sort {$a cmp $b} @{$tree_data{$tree_file}{$t_count}{$anc}{Target}{Species}};
			my @sc2a = sort {$a cmp $b} @{$tree_data{$tree_file}{$t_count}{$anc}{Other}{Species}};

			if($tree_data{$tree_file}{$t_count}{$anc}{Paralog} eq "Yes"){
					my @temp_sc = @sc1a;
					push(@temp_sc, @sc2a);
					@temp_sc = sort {$a cmp $b} @temp_sc;
					my %temp_ids;
					for my $sp (@temp_sc){
							$sp =~ /^(.*?)-\d/;
							$temp_ids{$1}=1;

					}
					my @names_temp;
					for my $k (keys %temp_ids){
							push (@names_temp, $k);
					}
					@names_temp = sort @names_temp;
					my $tsp = join (",", @names_temp);
					$final_counts{$tsp}++;
			}
			else{
					my @temp_sc1 = @sc1a;
					my @temp_sc2 = @sc2a;

					@temp_sc1 = sort {$a cmp $b} @temp_sc1;
					@temp_sc2 = sort {$a cmp $b} @temp_sc2;

					my %temp_ids1;
					my %temp_ids2;
					for my $sp (@temp_sc1){
							$sp =~ /^(.*?)-/;
							$temp_ids1{$1}=1;

					}
					for my $sp (@temp_sc2){
							$sp =~ /^(.*?)-/;
							$temp_ids2{$1}=1;

					}

					my @names_temp1;
					my @names_temp2;
					for my $k (keys %temp_ids1){
							push (@names_temp1, $k);
					}
					for my $k (keys %temp_ids2){
							push (@names_temp2, $k);
					}
					@names_temp1 = sort @names_temp1;
					@names_temp2 = sort @names_temp2;
					my $tsp1 = join (",", @names_temp1);
					$final_counts{$tsp1}++;

					my $tsp2 = join (",", @names_temp2);
					$final_counts{$tsp2}++;
			}


			

		my $sc1 = join (",", @{$tree_data{$tree_file}{$t_count}{$anc}{Target}{Species}});
		my $sc2 = join (",", @{$tree_data{$tree_file}{$t_count}{$anc}{Other}{Species}});

		print $out "$tree_file\t$tree_data{$tree_file}{$t_count}{$anc}{BS}\t$tree_data{$tree_file}{$t_count}{$anc}{BOTH}\t$tree_data{$tree_file}{$t_count}{$anc}{Target}{BS}\t$tree_data{$tree_file}{$t_count}{$anc}{Paralog}\t$sc1\t$tree_data{$tree_file}{$t_count}{$anc}{Other}{BS}\t$sc2\n";
		}}
	}
	for my $tfname (keys %missing_pars){
		$count_tree+=$missing_pars{$tfname};
	}

	my $split_partree=0;
	for my $sfname (keys %split_pars){
		$split_partree+=$split_pars{$sfname};
	}
	open my $out_f, ">", "Subgenome_Counts_" . $name . ".txt";
	for my $kid (sort {$a cmp $b} keys %final_counts){
			print $out_f "$kid\t$final_counts{$kid}\n";
	}
	print "Total paralogs found: $count_par\nTrees with paralogs: $count_tree\nTrees with Split Paralogs:$split_partree\n";
}
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

sub get_descendents_boot {
	#takes array of descendents and ancestor
	my $array_ref = $_[0];
	my @old_descendents= $array_ref;
	my $ancestor = $_[1]->ancestor;
	#my $bs = $ancestor->bootstrap;
	#if(!$bs){
	#		$bs = 0;
	#}
	my @temp_descendents = $ancestor->get_all_Descendents;
	my @descendents=();
	for my $temp_descend (@temp_descendents){
		if($temp_descend->is_Leaf){
			$temp_descend=$temp_descend->id;
			push(@descendents, $temp_descend);
		}
	}	
	return (\@descendents, $ancestor);

}
