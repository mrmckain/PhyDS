#!/usr/bin/perl
use warnings;
use strict;

my %counts;
open my $file, "<", $ARGV[0];
while(<$file>){
	chomp;
	my $spid;
	if($_ =~ /Panicum|Setar|Dichan|Andro|Sbi|Arun|Misc|Sacc|Them|Both|Schiz/){
		next;
	}
	if($_ =~ /,/){
		my @sparray = split ( /,/,$_);
		my $zea=0;
		my $uv=0;
		my $hemco=0;
		my $chio=0;
		my $thelo=0;
		my $arth=0;
		my $chryso=0;
		for my $sp(@sparray){
			if($sp =~ /Zmays|Zea|Tripsacum/){
				$zea++;
			}
			if($sp =~ /Vossia|Urel/){
				$uv++;
			}
			if($sp =~ /Hemar|Coel/){
				$hemco++;
			}
			if($sp =~ /Chio/){
				$chio++;
			}
			if($sp =~ /Thelo/){
				$thelo++;
			}
			if($sp =~ /Arthra/){
				$arth++;
			}
			if($sp =~ /Chryso/){
				$chryso++;
			}
		}
		if($zea == scalar @sparray){
			$spid = "Tripsacum,Zea,Zmays";
			$counts{$spid}++;
		}
		
		if($uv > 0 && $chio == 0 && $thelo == 0 && $hemco == 0 && $arth ==0){
			$counts{"Urelytrum,Vossia"}++;
		}
		
		if($uv == 0 && $chio > 0 && $thelo == 0 && $hemco == 0 && $arth ==0){
			$counts{"Chionachne"}++;
		}
		
		if($uv == 0 && $chio == 0 && $thelo > 0 && $hemco == 0 && $arth ==0){
			$counts{"Thelopogon"}++;
		}
		
		if($uv == 0 && $chio == 0 && $thelo == 0 && $hemco > 0 && $arth ==0){
			$counts{"Coelorachis,Hemarthria"}++;
		}
		
		if($uv == 0 && $chio == 0 && $thelo == 0 && $hemco == 0 && $arth >0){
			$counts{"Arthraxon"}++;
		}
		if($uv == 0 && $chio == 0 && $thelo == 0 && $hemco == 0 && $arth == 0 && $chryso >0){
			$counts{"Chrysopogon"}++;
		}
	
		}
		else{
		my $sp = $_;
		my $zea=0;
		my $uv=0;
		my $hemco=0;
		my $chio=0;
		my $thelo=0;
		my $arth=0;
		my $chryso=0;
		
			if($sp =~ /Zmays|Zea|Tripsacum/){
				$zea++;
			}
			if($sp =~ /Vossia|Urel/){
				$uv++;
			}
			if($sp =~ /Hemar|Coel/){
				$hemco++;
			}
			if($sp =~ /Chio/){
				$chio++;
			}
			if($sp =~ /Thelo/){
				$thelo++;
			}
			if($sp =~ /Arthra/){
				$arth++;
			}
			if($sp =~ /Chryso/){
				$chryso++;
			}
		
		if($zea > 0){
			$spid = "Tripsacum,Zea,Zmays";
			$counts{$spid}++;
		}
		
		if($uv > 0 && $chio == 0 && $thelo == 0 && $hemco == 0 && $arth ==0){
			$counts{"Urelytrum,Vossia"}++;
		}
		if($uv == 0 && $chio > 0 && $thelo == 0 && $hemco == 0 && $arth ==0){
			$counts{"Chionachne"}++;
		}
		
		if($uv == 0 && $chio == 0 && $thelo > 0 && $hemco == 0 && $arth ==0){
			$counts{"Thelopogon"}++;
		}
		
		if($uv == 0 && $chio == 0 && $thelo == 0 && $hemco > 0 && $arth ==0){
			$counts{"Coelorachis,Hemarthria"}++;
		}
		
		if($uv == 0 && $chio == 0 && $thelo == 0 && $hemco == 0 && $arth >0){
			$counts{"Arthraxon"}++;
		}
		if($uv == 0 && $chio == 0 && $thelo == 0 && $hemco == 0 && $arth == 0 && $chryso >0){
			$counts{"Chrysopogon"}++;
		}
		
		}
	}

open my $out, ">", "Count_Zmm_HIT_Outside_Taxa.txt";

for my $sp (keys %counts){
	print $out "$sp\t$counts{$sp}\n";
}
	
	