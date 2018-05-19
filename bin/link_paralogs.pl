#!/usr/bin/perl -w
use strict;

my %trees;
open my $treefile, "<", $ARGV[2];
while(<$treefile>){
	chomp;
	my @tarray = split /\s+/;
	my @gene_id;
	if($tarray[5] =~ /,/){
		my @g_array = split(/,/, $tarray[5]);
		for my $gzt (@g_array){
			if($gzt =~ /Zmays.(.*?)$/){
				push (@gene_id, $1);
			}
		}
		#@gene_id = ($tarray[5] =~ /Zmays.(.*?),Zmays.(.*?)$/g);
		#$gene_id=$1;
	}
	else{
		
		$tarray[5] =~ /Zmays.(.*?)$/;
		push(@gene_id, $1);
		#$gene_id = $1;
	}
	
	$tarray[0] =~ /Urel_voss\/RAxML_bipartitions.(.*?).cds_cluster..fasta.raxml.out/;
	my $ortho = $1;
	for my $val (@gene_id){
		#$val =~ /Zmays.(.*?)$/;
		$trees{$val}=$ortho;
	}
}

my %paralogs;
open my $paralogfile, "<", $ARGV[0];
while(<$paralogfile>){
        chomp;
        my @tarray= split /\s+/;
        $tarray[0] =~ /Zmays.(.*?)$/;
        $tarray[0]=$1;
        $tarray[1] =~ /Zmays.(.*?)$/;
        $tarray[1]=$1;
        $paralogs{$tarray[0]}=$tarray[1];
        $paralogs{$tarray[1]}=$tarray[0];
}

my %hits;
open my $file, "<", $ARGV[1];
while(<$file>){
	chomp;
	if(/Chromosome/){
		next;
	}
	my @tarray=split/\s+/;
	$hits{$tarray[1]}=$_;
}
		
my %used;
open my $out, ">", "Linked_subgenomes_via_syntelogs_BSV$ARGV[3].txt";
print $out "Chromosome\tGene_ID\tStart\tPutative_Sub\tOrthogroup\tChromosome\tGene_ID\tStart\tPutative_Sub\tOrthogroup\n";
for my $gene (keys %hits){
	if(exists $used{$gene}){
		next;
	}
	print "$gene\n";
	print "Tree= $trees{$gene}\n";
	if(exists $paralogs{$gene} && exists $hits{$paralogs{$gene}}){
		my @tarray1 = split(/\s+/, $hits{$gene});
		my @tarray2 = split(/\s+/, $hits{$paralogs{$gene}});
		if($tarray1[0] < $tarray2[0]){
			print $out "$hits{$gene}\t$trees{$gene}\t$hits{$paralogs{$gene}}\t$trees{$paralogs{$gene}}\n";
		}
		else{
			print $out "$hits{$paralogs{$gene}}\t$trees{$paralogs{$gene}}\t$hits{$gene}\t$trees{$gene}\n";
		}
			
		$used{$gene}=1;
		$used{$paralogs{$gene}}=1;
	}
	else{
		print $out "$hits{$gene}\t$trees{$gene}\tNA\tNA\tNA\tNA\tNA\n";
}
}


