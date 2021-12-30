#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::TreeIO;
use Getopt::Long;
use Pod::Usage;

my $debug = 0;
my $man = 0;
my $help = 0;
my $paralogs;
my $trees;
my $outgroups;
my $species_tree;
my $prefix = "PhyDS";
my $all_genes;
my $estimate_paralogs;
my $tree_type="ML";

GetOptions('help|?' => \$help, man => \$man, "trees=s" => \$trees, "outgroups=s" => \$outgroups, "name=s" => \$prefix) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


=head1 NAME
PUG -- Phylogentic Placement of Polyploidy Using Genomes
=head1 SYNOPSIS
perl gapid.pl -trees file -outgroups list [options] 
=head1 OPTIONS
    -trees             Directory of gene trees used to identify WGD placement. Trees should be in Newick format and have bootstrap values. 
    -outgroups         Comma delimited list of outgroups for trees. At least one of these taxa must be in the tree for it to be considered.
    -name              Identifier for this run.  [Default = "PhyDS"]
    -man               Full documentation
    -debug             Turn on debugging messages
=cut


my %results;
my %dups;
my %dups2;
my %seqpairs;
my %seqpairs2;
my %lcas;
my %events;
my %hypotheses;
my %species_taxa;
my %species_index;
my %node_index;
my $node_count=0;


###Estimate all possible pairs if option chosen.###
my @file = <$trees/*>;
my %putative_paralogs;
open my $out_estparalogs, ">", "$prefix\_Estimated_Putative_Paralogs.txt";
$paralogs = "$prefix\_Estimated_Putative_Paralogs.txt";

for my $treefile (@file){
    my $treeio = new Bio::TreeIO(-format => "newick", -file => "$treefile", -internal_node_id => 'bootstrap');
    while( my $tree = $treeio->next_tree ) {
        my %temp_paralogs;
        my @taxa = $tree->get_leaf_nodes;
        for my $taxa (@taxa){
            my $taxon = $taxa->id;
            $taxon =~ /(.*?)-(.*?)/;
            $species_taxa{$1}=1;
        }
    }
}

for my $treefile (@file){
    my $treeio = new Bio::TreeIO(-format => "newick", -file => "$treefile", -internal_node_id => 'bootstrap');
    while( my $tree = $treeio->next_tree ) {
        my %temp_paralogs;
        my @taxa = $tree->get_leaf_nodes;
        foreach my $taxa (@taxa){
              my $taxon = $taxa->id;
             for my $spec_tax (keys %species_taxa){
                     if($taxon =~ /$spec_tax/){
                            $temp_paralogs{$spec_tax}{$taxon}=1;
                       }
             }
         }

        for my $ospecies (keys %temp_paralogs){
            if(scalar keys %{$temp_paralogs{$ospecies}} >= 2){
                    my @transcripts = keys %{$temp_paralogs{$ospecies}};
                    for (my $i=0; $i < scalar @transcripts-1; $i++){
                        for (my $k = $i+1; $k <= scalar @transcripts-1; $k++){
                            print $out_estparalogs "$transcripts[$i]\t$transcripts[$k]\tunknown\n";
                         }
                    }
            }
        }
    }
}

