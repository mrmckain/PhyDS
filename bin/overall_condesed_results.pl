#!/usr/bin/perl
use strict;
use warnings;

my @files = <Results_for_BSV*_unlinked_condensed*.txt>;

my %results;

for my $tfile (@files){
        open my $file, "<", $tfile;
        $tfile =~ /BSV(\d.*?)_unlinked_condensed/;
        my $bsv = $1;
        while(<$file>){
                chomp;
                my @tarray=split/\s+/;
                $results{$tarray[0]}{$bsv}=$tarray[1];
        }
}


open my $out, ">", "SG_ID_All_identity.txt";
print $out "Taxon\tBSV5\tBSV10\tBSV15\tBSV20\tBSV25\tBSV30\tBSV35\tBSV40\tBSV45\tBSV50\tBSV55\tBSV60\tBSV65\tBSV70\tBSV75\tBSV80\tBSV85\tBSV90\tBSV95\tBSV100\n";
for my $sp (sort {$a cmp $b} keys %results){
        print $out "$sp";
        for my $bb (sort {$a <=> $b} keys %{$results{$sp}}){
                
                print $out "\t$results{$sp}{$bb}";
        }
        print $out "\n";
}