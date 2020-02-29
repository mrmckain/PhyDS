PhyDS: Phylogenetic iDentification of Subgenomes 
=============
<b>Author</b>: Michael R. McKain<br>
</br>
Version 2.1<br>
</br>
<b>Contact</b>: https://github.com/mrmckain

<h1>Description</h1>
PhyDS is a series of scripts designed to identify the evolutionary origin of subgenomes. 

PhyDS takes a focal set of paralogs derived from either a genome (e.g. syntenic homoeologs) or transcriptome and searches a set of gene trees to identify what the most closely related gene sequences are to these focal paralogs. PhyDS first identifies the paralogs present in a gene tree and then moves to the direct ancestral node of the paralog. PhyDS then identifies all gene sequences (tips) that are direct descendants of the current node. Taxon identities of these sequences are assessed to determine if one of two conditions is met: 1) taxa are present that are not in a user-given list of taxa to ignore and 2) the paralog mate of the current focal paralog is found. If neither of these conditions is met, PhyDS will proceed to the ancestor node of the current node and repeat the process. Once these conditions are met, PhyDS records the bootstrap value of the node, the paralogs of interest, and the sequence identities found in the clade descendant from the current node.

<h1>Dependencies</h1>

PhyDS is a Perl-based program and was written using Perl 5.x. Two Perl modules are needed: Bio::TreeIO and Getopt::Long. Bio:TreeIO is a part of BioPerl and may need to be installed in your system. Getopt::Long is a standard module is likely already available.

<h1>Installation</h1>

No installation is required other than the above Perl modules.

<h1>Input</h1>

<h3>Paralogs</h3>

A tab- or space-delimited list of putative paralogs from your target taxon is needed for searching trees with PhyDS. Names of the sequences must match those used in the gene trees.

<h3>Ignore Taxa</h3>

This is a comma-separated list of taxa that PhyDS can ignore when identifying a stopping point to assigned the closest relative to a paralog. If you know that a few of your species share the same polyploidy event, then you would use a list of unique identifiers used in your tip (or sequence) names to tell PhyDS that if it finds these taxa that they are to be ignored. PhyDS is looking for either the putative paralog to your target paralog in the search or it is looking for a taxon that is not known to share the polyploidy event. When it finds either of these, it stops searching and returns the taxa mostly closely related to your target paralog. 

Example: If using <i>Zea mays</i> as a target taxon but also have <i>Zea_luxurians</i>, <i>Zea_perennis</i>, <i>Tripsacum_dactyloides</i>, <i>Tripsacum_peruvianum</i>, and <i>Tripsacum australe</i> in your gene trees, you can either use "Zea,Tripsacum" or "Zea_luxurians,Zea_perennis,Tripsacum_dactyloides,Tripsacum_peruvianum,Tripsacum_australe" for this option since the polyploidy event of maize is shared by both <i>Zea</i> and <i>Tripsacum</i>. You could also use "luxurians,perennis,dactyloids,peruvianum,australe" since they are unique aspects of the taxa names. These identifiers can be anywhere in the name of the sequence/tree tip.

<h3>Gene Trees</h3>

Gene trees are expected to be the most likely trees with bootstrap values. For example, the "bipartitions" file produced by RAxML works here. Other programs should be able to provide similarly formatted gene trees. All trees must be in the Newick format.

<h1>Output</h1>

To be completed.

<h3>General Syntax</h3>

<code>perl phyds.2.pl -p <paralogs_file> --trees <directory_of_gene_tree> --ignore <comma_separated_list_of_taxa> --name <sample_name> --bootstrap <minimum_bootstrap_cutoff> </code>

