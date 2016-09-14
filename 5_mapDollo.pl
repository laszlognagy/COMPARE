#!/usr/bin/perl

#	USAGE ./script --in <input nexus file> --spectree <organismal phylogeny file>
#	Will print to STDOUT
#	Note the $regex variable need to be edited on line 21 to conform to current naming scheme of proteins

use strict;
use warnings;
use Bio::Phylo::IO; 
use Bio::Phylo::Matrices;
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::Factory;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::TreeIO; #onnly eats rooted trees
use Getopt::Long;

my ($CharNo, $infile, $SpeciesTreeFile);
my (%TermNames, %Gains, %Losses, %Comments, %Matrix, %AncStates, %Cluster_lookup);
my @Names_in_subtree = my @nodestovisit = ();
my $regex = '(_\w+)'; #global variable accessible everywhere, no need to supply to subs

###############Parsing command line arguments############
GetOptions( 'in=s' => \$infile, #Working Directory
			'spectree=s' => \$SpeciesTreeFile); # 
#			'min:i' => \$MinSize, #minimum size of cluster to subject to blast-based check
#			'merged:s' => \$allgenomes, #merged genomes fasta file
#			'phyd:s' => \$phylodir, #single copy directoiry
#			'restart' => \$restart); # whether to restart or not

# Check command line parameters
if (not $SpeciesTreeFile or not (-e $SpeciesTreeFile)) { die "ERROR: missing input organismal tree file\n"; }
if (not $infile or not (-e $infile)) { die "ERROR: missing input nexus file\n"; }
sleep(1);

#when correlating with trait presence, the both gains and losses for a single character have to be represented in one single data structure

################# Parsing species tree and input NEXUS file ##################
my $SpeciesTree = Bio::Phylo::IO->parse(
   -file => "$SpeciesTreeFile",
   -format => 'newick'
)->first;				#probably we should take care to involve phylogenetic uncertainty here by iterating the whole process over several species trees?????????

my $No_of_Chars = &parse_nexus($infile);

#foreach (keys %Matrix) { if you want the full matrix in a more readable form, uncomment this
#	print "matrix $_ @{$Matrix{$_}}\n";
#}

foreach my $terminal (my @Terminals = @{$SpeciesTree->get_terminals}) {
	my $id = $terminal->id;
	my $intid = $terminal->internal_id;
	$TermNames{$id} = $terminal; #taxonnames(human readable) x internal Node object
}
##############


#############MAIN#############
for my $i (0 .. $No_of_Chars - 1) { #should start at 1 since the character numbering starts at one oalso ($#{$Matrix{int rand keys %Matrix}})
	print "processing character $i\n";
	%AncStates = ();
	my %SingleChar;
	foreach (keys %Matrix) { #put together the singlechar hash
		$SingleChar{$_} = $Matrix{$_}[$i];
		$CharNo = $i + 1;
	}
	###########Finding the place of gain on the tree##########
	my @Terminals_in_subtree = (); #ez kiemelheto innen valahova
	foreach (@{$Comments{$i+1}}) {
		push(@Terminals_in_subtree, $TermNames{$_});
	}
	my $mrca = $SpeciesTree->get_mrca(\@Terminals_in_subtree);
			my $mrcaid = $mrca->internal_id;
			my @mrcadesc = @{$mrca->get_terminals};
	
	my $mrca_id = $mrca->internal_id;
	my $info = "$Cluster_lookup{$CharNo}" . '/' . "$CharNo";
	push (@{$Gains{$mrca_id}}, $info);
	#########################


	############Reconstruct loss(es)################

	my @Descendents = $mrca->get_all_Descendents;
	push (@Descendents, $mrca);
	foreach (@Descendents) {
		if ($_->is_terminal) {
			my $id = $_->id;
			$AncStates{$_} = $SingleChar{$id}; #Ancstate is node_object x character state. for internal nodes no state is assigned at this point
		} else {
			$AncStates{$_} = ();
		}
	}
	@nodestovisit = @Descendents; #this line is superfluous
	until ((scalar @nodestovisit) == 0) {
		my $node = shift @nodestovisit;
		next if ($node->is_terminal);
		my $returnakarmi = &traverser(\$node);
	}
	foreach (@Descendents) { #this foreach is only for checking purposes
		next if ($_->is_terminal);
		my $id = $_->internal_id;
		my (@termlist, @Ters);
		@Ters = @{$_->get_terminals};
		foreach my $te (@Ters) {
			my $ud = $te->id;
			push (@termlist, $ud);
		}
	}
	
	##############Find state changes along the tree######
	$mrca->visit_level_order(sub{ 
				my $node = shift;
				my $parent = $node->get_parent;
				my $rootnode = $SpeciesTree->get_root();
				if ($node ne $rootnode) {
					my $try = $parent->internal_id;
					if (exists $AncStates{$parent}) { #with this we rule out the parent of $mrca ########oke de mi van a rootnode-nal?
						if ($AncStates{$node} ne $AncStates{$parent}) { #in this case there was a loss along that node
							my $intid = $node->internal_id; #terminalisokra ez mukodik?
							my $info = "$Cluster_lookup{$CharNo}" . '/' . "$CharNo";
							push (@{$Losses{$intid}}, $info); #this way, for each node on the tree, we will have one array in the hash which contains all the character numbers which show a loss
						}
					}
				}
	});
}	
############End of Main############

#############print results
#print "final szamvetes\n";

my $rootnode = $SpeciesTree->get_root(); #ez a resz ebben a scriptben nem ,lett kiprobalva, de attol mukodnie kellene schiman
print "For your reference, here is the list of internal nodes of the tree and the descendants for each\n";
$rootnode->visit_level_order(sub{
	my $node = shift;
	my @terminals = ();
	my $nodename = $node->internal_id;
	foreach (my @nodes = @{$node->get_terminals}) {
		my $terID = $_->id;
		push(@terminals, $terID);
	}
	print "node$nodename\t @terminals\n";

});



foreach (sort {$a <=> $b } keys %Gains) {
	print "GAIN\tnode\t$_\tcharacters\t@{$Gains{$_}}\n";
}
foreach (sort {$a <=> $b } keys %Losses) {
	print "LOSSES\tnode\t$_\tcharacters\t@{$Losses{$_}}\n";
}


#############SUBS###############

sub traverser {
	my $startnode_r = shift;
	my $startnode = $$startnode_r;
	my @kids = @{$startnode->get_children};
	if (defined $AncStates{$kids[0]} and defined $AncStates{$kids[1]}) { #if there is state info for both immediate chldren
		if (!defined $AncStates{$startnode}) { #ez ujide be kell tenni egyet ami kivedi hogy mar rekonstrualt noduszok ujra rekonstrualodjanak next if defined ancsates{f}
			&ancstate($startnode);#lehet gyorsabb lenne betolni neki a @kids arrayt mint ujra kiszedni
		}
	}else {
		push (@nodestovisit, $startnode); #here the node is placed back into the nodestovisit array, because the ancstate could not be determined
		foreach (@kids) {
			if ($_->is_internal) {
				next if (defined $AncStates{$_});
				&traverser(\$_);
			}
		}
	}
}


sub ancstate {
	my $startnode = shift;
	my $nodeid = $startnode->internal_id;
	my @kids = @{$startnode->get_children};
	if (($AncStates{$kids[0]} + $AncStates{$kids[1]}) >= 1) { #in this case the ancstate must be 1
		$AncStates{$startnode} = '1';
	} else {
		$AncStates{$startnode} = '0';
	}
}


sub regex {
	my $in = shift;
	my $tID = ();
	my @array = @$in;
	my @ofspecnames = ();
	foreach (@array) {
		if ($_ =~ /$regex/i) {
				$tID = $1;
				$tID =~ s/_//;
		}
		push(@ofspecnames, $tID);
		}
	my @uns = keys %{{ map{$_=>1}@ofspecnames}};
	return(\@uns);
}

sub parse_nexus {
	my $inf = shift;
	my @Names = ();
	my $i = 1;
	my $No_of_Chars = ();
	my @Terminals = @{$SpeciesTree->get_terminals};
	foreach (@Terminals) {
		my $name = $_->id;
		push(@Names, $name);
	}
	my %hash = map { $_ => 1 } @Names;
	open INFILE, "<$inf" || die "infile could not be opened: $!\n";
	while (<INFILE>) {
		chomp;
		if ($_ =~ /NCHAR=(\d*);/) {
			$No_of_Chars = $1;
		}
		if ($_ =~ /^\s?(\w*)\s{10}/) { #Storing matrix in %Matrix HoA
			if (exists $hash{$1}) {
				$1 =~ s/\s?$1\s{10}//; #strip out 
				my @Line = split(' ',$_);
				shift @Line;
				push @{$Matrix{$1}}, @Line;
			}	
		}
		if ($_ =~ /\[\d/) { #Storing comments (proteins and sistertaxa info) in %Comments
			s/\]$//;
			$_ =~ /Cluster(\d*).*protein\(s\)(.*)\t\t\t\tSistertaxa/; #[317: is coded from gene tree ./Tree_files_rerooted/Cluster4.fastaaligned.out.2.fas.tree.out.tre protein(s) 50673|Dacsp1 2413|Sporo1 346907|Serla_S7_9 143804|Conpu1 17021|Wolco1 150347|Dicsq1 55517|Stehi1 3147|Malgl1 139579|Phchr1 142066|Fompi1 118145|Trave1 230657|Schco1 48594|Hetan2 113950|Fomme1 136215|Glotr1_1 6895|Cryne_H99_1 7465|Stano1				Sistertaxa: Punst1 Picst3 Mellp1 Aurde1 Trire2 Ustma1 Aspni5 Malgl1 Stano1]
			my $ClustNo = $1;
			my @Prots = split(' ', $2); #proteins that make up the character
			my $Prots_species_ref = &regex(\@Prots);
			my @Prots_species = @$Prots_species_ref;
			$_ =~ /Sistertaxa:\s(.*)/;
			my @Sis = split(' ', $1); #with the sister taxa will the set of species that has the gene be ready
			push(@Sis, @Prots_species);
			@Names_in_subtree = keys %{{ map{$_=>1}@Sis}};#eliminate duplicates
			@{$Comments{$i}} = @Names_in_subtree; #assigning commentstring to charnumber
			$Cluster_lookup{$i} = $ClustNo;
			$i++;
		}
	}
	close INFILE;
	return ($No_of_Chars);
}
