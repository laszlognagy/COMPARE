#!/usr/bin/perl

#	USAGE ./script --tree_d <path to trees(e.g. /users/bob/project/trees or ./Treefiles)> --spectree <organismal phylogeny file> --clusterfile <clusterfile> [--exclude <file containing clusters to exclude>] [--outfile <outfilename>]
#	Note the $regex variable need to be edited on line 20 to conform to naming scheme of proteins
#	In the current implementation, the script looks for reconciled tree files with the extension '.tree' (treefix default). Trees should be a single line (unlike treefix default). Files are assumed to start with Cluster<number>....


#use strict;
use warnings;
use Bio::Phylo::IO; 
use Bio::Phylo::Matrices;
use Bio::Phylo::Matrices::Matrix;
use Bio::Phylo::Factory;
use Bio::Tree::Tree;
use Bio::Tree::TreeFunctionsI;
use Bio::TreeIO; #onnly eats rooted trees
use Getopt::Long;

my $i = 0; 
my $regex = '(_\w+)';
my ($decision, $taxonID, $index);
my (@SpecIDs, @DuplNodes, @SpecNodes, @Proteins_done, @Nodes_to_Visit, @allnodes);
my (%Matrix, %Comments, %FirstPairs, %Sisters, %LastPairs, %goodnodes, %Exclusion);


###############Parsing command line arguments############
GetOptions( 'tree_d=s' => \$tree_dir, #Working Directory
			'rt=s' => \$rooting_method, 
			'spectree=s' => \$SpecTreeFile, #
			'outfile=s' => \$OutFileName, #
			'exclude:s' => \$BlastFileName, #
			'clusterfile=s' => \$ClustFile); #Clustering file

# Check command line parameters
if (not $SpecTreeFile or not (-e $SpecTreeFile)) { die "ERROR: missing input organismal tree file\n"; }
if (not $tree_dir){ 
	$tree_dir = ".";
	print "WARNING: You have not specified a tree directory, so I assume your input newick files are in the current folder\n";
}
if (not $ClustFile){ 
	$ClustFile = "";
	print "WARNING: You have not specified a cluster file, so I assume you dont want to recode cluasters smaller than 4 proteins\n";
}
if (not $OutFileName){ 
	print "WARNING: You have not specified any output file, results will be printed to Matrix.nex\n"; 
	$OutFileName = 'Matrix.nex';
}
#if (not $rooting_method){ 
#	print "You have not specified how to root the trees, midpoint rooting will be used\n"; 
#}
print "$SpecTreeFile $tree_dir\n";

if ($BlastFileName) {
	open EXCLUDE, "<$BlastFileName" || die "the file you specified for excluding clusters could not be found\n";
	while (<EXCLUDE>) {
		if ($_ =~ /Cluster(\d*)/) {
			$Exclusion{$1} = $_;
		}
	}
}

sleep(1);

####a feherjenevek regexes kiszedes 9x van meg, ki lehetne tenni sub-ba

################# Parsing species tree ##################
my $SpeciesTree = Bio::Phylo::IO->parse(
   -file => "$SpecTreeFile",
   -format => 'newick'
)->first;				#first says first tree in the file

foreach my $a (my @SpecLabels = @{$SpeciesTree->get_terminals}) {
	my $id = $a->id;
	push(@SpecIDs, $id);
}
my $nntax = scalar @SpecIDs;
##########################

################# Main ##################
print "starting the magic in directory $tree_dir\n";
my @TreeList = <$tree_dir/*>;
foreach (@TreeList) {
	if ($_ =~ /Cluster(\d*)\./) {
		$Trees{$1} = $_;
	}
}

	################# Parsing gene trees ##################
foreach (my @SortedTrees = sort {$a <=> $b } keys %Trees) { #for each gene tree in folder
	my $single = $Trees{$_};
	if ($single =~ /.tree/) {
		next if (defined $Exclusion{$_});
		print "Processing $single\n";
		undef @Matrix{@SpecIDs};	# set an undef value for each species of the species tree in the big matrix
		undef %Comments;
		undef %Sisters;
		@Nodes_to_Visit = @allnodes = @Proteins_done = ();
		#here %Matrix and Comment and third should be emptied
		my $tree = Bio::Phylo::IO->parse(
		  -file => "$single", 
		  -format => 'newick'
		)->first;
		my $rootnode = $tree->get_root;
		push (@Nodes_to_Visit, $rootnode);
		push (@allnodes, $rootnode);
		foreach (@Nodes_to_Visit) {
			@SpecNodes = @DuplNodes = ();
			my @backtrack = ();
			undef %SpecNodes_h;
			undef %FirstPairs;
			$stolennode = shift @Nodes_to_Visit;
			push (@allnodes,$stolennode);
			&FirstScan($stolennode); # look for the path containing the most speciation events across the tree
			if ($#SpecNodes >= 0)  { #uj
				foreach (@SpecNodes) {
					my $u = $_->internal_id;
					$SpecNodes_h{$u} = $_;
				}
					@backtrack = reverse sort {$a <=> $b} keys %SpecNodes_h;
					my $ancestry_ref = &ancestry(\@backtrack, \$tree); 
					my @ancestry = @$ancestry_ref;
					$firstpairs_ref = &backtracker(\@backtrack);
					%FirstPairs = %$firstpairs_ref;
					if ((scalar keys %FirstPairs) > 0) {
						&code_it(\%FirstPairs, \@ancestry);
					}
} 
		}


		$tree->visit_level_order(sub{ 
			my $node = shift;
			if (!grep {$_ eq $node} @allnodes) {
				my $nnnode = $node->internal_id;
				if ($node->is_internal) {
						push(@allnodes, $node);
						@backtrack = @SpecNodes = @DuplNodes = ();
						undef %SpecNodes_h;
						undef %FirstPairs;
						&FirstScan($node);
						foreach (@SpecNodes) {
							my $u = $_->internal_id;
							$SpecNodes_h{$u} = $_;
						}
						@backtrack = reverse sort {$a <=> $b} keys %SpecNodes_h;
						$firstpairs_ref = &backtracker(\@backtrack);
						%FirstPairs = %$firstpairs_ref;
						my $ancestry_ref = &ancestry(\@backtrack, \$tree);
						my @ancestry = @$ancestry_ref;
						if ((scalar keys %FirstPairs) > 0) {
							&code_it(\%FirstPairs, \@ancestry);
						}
				}
			}
		});

		############ Final traversal ###########
		foreach $f (my @allterminals = @{$rootnode->get_terminals}) {
			undef %LastPairs;
			$lastid = $f->id;
			if (!grep {$_ eq $lastid} @Proteins_done) {
				if ($lastid =~ /$regex/i){  #28580|Mellp1
					$lID = $1;
					$lID =~ s/_//; 
				}
				push (@Proteins_done, $lastid);
				$LastPairs{$lID} = $lastid;
				my $Sistertaxa_r = &last_ancestry(\$f);
				my @Sistertaxa = @$Sistertaxa_r;
				&code_it(\%LastPairs, \@Sistertaxa);
			}
		}
		($i, $nchar) = &Dumper($single, $OutFileName, $i); #, \%Matrix, \%Comments
	}
}
		$nchar = &small_clusters($ClustFile, \%Exclusion);
		&finish_Matrix($OutFileName, $nchar); # adding NCHAR and endblock

########################### end of big foreach 

#############################Subs################################
sub FirstScan { #this one runs along the tree to find the path that maximizes the number of speciation events. Proteins along this will be designated as the survivors of the ancestral paralog
	my $curr_node = shift;
	my @terminalslist = @{$curr_node->get_terminals}; #kell ez ide egyaltalan????
	foreach my $elem (@terminalslist) {
		my $tu = $elem->id;
	}
	my $winner;
	if ($curr_node->is_preterminal) {
		push(@SpecNodes, $curr_node);
	} else {
		my @childs = @{$curr_node->get_children};
		foreach my $childs (@childs) {
			my $nev = $childs->internal_id;
			my $nevv = $childs->id;
			next if ($childs->is_terminal); #ez biztos jo, le lett tesztelve
			if (!grep {$_ eq $childs} @allnodes) {
				my ($decision) = &compare_subtrees($curr_node); #check subtree by sub
				if ($decision == 0) { #speciation node		
					$winner = $childs;
					if ($winner->is_internal) {
						if (!@SpecNodes || !grep  {$_ eq $curr_node} @SpecNodes) {
							push(@SpecNodes, $curr_node);
						}
					}
				}
				if ($decision == 1) { #if duplication node
					$winner = &which_subtree_has_more_species(\$curr_node);
					my $nev = $winner->internal_id;
			my $nevv = $winner->id;
			my @terminallist = @{$winner->get_terminals};
	foreach my $velem (@terminallist) {
		my $tru = $velem->id;
	}
					if ($winner->is_internal) {
						if (!@DuplNodes || !grep  {$_ eq $winner} @DuplNodes) {
							push(@DuplNodes, $curr_node);
						}
						if (!@DuplNodes || !grep  {$_ eq $curr_node} @DuplNodes) {
							push(@DuplNodes, $curr_node);
						}
					}
					my $cladesize_r = &is_monotypic(\$curr_node);
					my $cladesize = $$cladesize_r;
					if ($cladesize == 1) {
						push(@SpecNodes,$curr_node);
						return;
					}
				}
				if ($winner->is_internal) {
					push(@allnodes,$winner);
					&FirstScan($winner);
				}
			}
		}
	}
}

sub which_subtree_has_more_species {
	my $inputnode_ref = shift;
	my $inputnode = $$inputnode_ref;
	my $tID;
	my %uniques;
	my @subchilds = $inputnode->each_Descendent();
	foreach $subchilds (@subchilds) {
		undef @ofspecnames;
		foreach $kids (my @kids = @{$subchilds->get_terminals}) {
			my $ezt = $kids->id;
			if ($ezt =~ /$regex/i){
				$tID = $1;
				$tID =~ s/_//;
			}
			push(@ofspecnames, $tID);
		}
		 @uns = keys %{{ map{$_=>1}@ofspecnames}};
		 $uniques{$subchilds} = scalar @uns;
	}
	my $subchilds0 = $subchilds[0];
	my $subchilds1 = $subchilds[1];
	if ($uniques{$subchilds0} > $uniques{$subchilds1}) {
		if (not $subchilds1->is_terminal) {
			push (@Nodes_to_Visit, $subchilds1);
		}
		return $subchilds0;
	}
	else {
		if (not $subchilds0->is_terminal) {
			push (@Nodes_to_Visit, $subchilds0);
		}
		return $subchilds1;
	}
}

sub backtracker {
	my $backtr_ref = shift;
	my @backtrack = @$backtr_ref;
	foreach $nodes (@backtrack) { # for each node along the path containing the most speciation nodes, get species in that subtree and add them to the given character in the matrix. Walk is from the leaves towards the root. 
		foreach $kolkok (my @kolkok = $SpecNodes_h{$nodes}->each_Descendent) { 
			if ($kolkok->is_terminal) {
				my $p = $kolkok->id; #$p is the prot is that should be stored
				if (!grep($_ eq $p, @Proteins_done)) {
					if ($p =~ /$regex/i){ 
						$dID = $1;
						$dID =~ s/_//; # this is for the format Speciesname|ProtID|whateverdetail
					}
					if (!grep  {$_ eq $dID} keys %FirstPairs) {
						push (@Proteins_done, $p);
						$FirstPairs{$dID} = $p;
					}
				}
			}
			if ($kolkok->is_internal) {
				if (!grep {$_ eq $kolkok} @SpecNodes) { #this ensures that only those childs of a speciation node will be visited here, which are not parts of the main backtrack path
					if (!grep {$_ eq $kolkok} @DuplNodes) { #this ensures that only those childs of a speciation node will be visited here, which are not parts of the main backtrack path
					my @SpecNodes_archive = @SpecNodes; #see 4 lines later
					my $difference = ();
					push(@allnodes, $kolkok);
					foreach my $SpN (@SpecNodes) { #this is to check if the &firstscan called above has changed anything the the specnodes array. If it has not changed anything, it means that both childs checked in firstscan were already present in the allnodes(?) array
						if (!grep {$_ eq $SpN} @SpecNodes_archive) {
							$difference++;
						}
					}
					foreach (@SpecNodes) {
						my $u = $_->internal_id;
						$SpecNodes_h{$u} = $_;
					}
					foreach $kids (my @kids = @{$SpecNodes_h{$nodes}->get_terminals}){ 
						my $p = $kids->id; #$p is the prot is that should be stored
						if  (!grep($_ eq $p, @Proteins_done)) {
						if ($p =~ /$regex/i){ 
							$kID = $1;
							$kID =~ s/_//; # this is for the format Speciesname|ProtID|whateverdetail
						}
						if (!grep  {$_ eq $kID} keys %FirstPairs) { #if species $kID is not yet present in the hash
							push (@Proteins_done, $p);
							$FirstPairs{$kID} = $p;
						}
						}
					}
					}
				}
			}
		}
	}
	return (\%FirstPairs);
}

sub compare_subtrees { #this takes descendants of two subtrees and checks if any of the species is shared between those subtrees
	my $current = shift;
	my (@desc_sis0, @desc_sis1, @tax_sis0, @tax_sis1);
	my @sisnodes = $current->each_Descendent(); #from here on we assume all trees are fully binary #this method crashes on polytomies?????
	@desc_sis0 = @{$sisnodes[0] ->get_terminals()}; #this might be solved easier with get_terminals
	@desc_sis1 = @{$sisnodes[1] ->get_terminals()};
	if ($sisnodes[0]->is_terminal) { 
		push(@desc_sis0, $sisnodes[0]);
	}
	if ($sisnodes[1]->is_terminal) {
		push(@desc_sis1, $sisnodes[1]);
	}
	foreach my $element (@desc_sis0) {
			my $tax_id = $element->id;
			if ($tax_id =~ /$regex/i){ 
				$taxonID = $1;
				$taxonID =~ s/_//; 
			}
			push(@tax_sis0, $taxonID); #here could come something that accounts for the naming scheme of proteins, which might be defined upon launching
	}
	foreach my $element (@desc_sis1) {
			my $tax_id = $element->id;
			if ($tax_id =~ /$regex/i){ 
				$taxonID = $1;
				$taxonID =~ s/_//; 
			}
			push(@tax_sis1, $taxonID);
	}
	my %temphash = map{$_ => 1} @tax_sis0;
	my $overlaps = scalar(grep($temphash{$_}, @tax_sis1));
	if ($overlaps > 0) {	
		$decision = 1; #this case its a duplication node
	}
	else {$decision = 0};
	return ($decision)
}


sub code_it {
	my $hashref = shift;
	my $ancref = shift;
	my %ProtPairs = %$hashref; 
	my @sistertaxa = @$ancref;
	my %SingleChar;
	my @Prot_IDs;
	foreach (keys %ProtPairs) {
		$SingleChar{$_} = 1; 
	}
	foreach (@SpecIDs) {
		next if (exists $SingleChar{$_});
		if (!exists $SingleChar{$_}) { 
			$SingleChar{$_} = 0; 
		}
	}
	foreach $key (keys %Matrix) {
		my @value = @Matrix{@key};
		push(@{$Matrix{$key}}, $SingleChar{$key});
	}
	foreach my $ids (values %ProtPairs) {
		push(@Prot_IDs, $ids);
	}
	push(@Proteins_done, @{values %ProtPairs});
	$index = (scalar keys(%Comments)) + 1; 
	my $CommentString = (join ' ',(values %ProtPairs));
	$Comments{$index} = $CommentString;
	my $sisterstring = (join ' ', @sistertaxa);
	$Sisters{$index} = $sisterstring;
}

sub Dumper {
	my $TreeFileName = shift;
	my $OutFile = shift;
	my $i = shift;
	my $tre;
	if ($i == 0) {
		if (-e $OutFile) { #removes already existing matrices with the same filename
			unlink $OutFile;
		}
	}
	my $newchar = scalar values %Comments;
	open OUTPUT, ">>$OutFile" || die "creation of output file failed: $!\n"; 
	if ($i == 0) {
		$nchar = 0;
		my $time = &mytime;
		my @taxa = keys %Matrix;
		print OUTPUT "#NEXUS\n[By convention every nexus file has to contain at least on line of comment: File written on $time by -Magnificent name to come- from all trees found in $tree_dir]\n\n\nBEGIN TAXA;\n\tTITLE Taxa;\n\tDIMENSIONS NTAX=$nntax;\n\tTAXLABELS\n\t\t@taxa\n\t;\n\nEND;\n\nBEGIN CHARACTERS;\n\tTITLE  Character_Matrix;\n\tDIMENSIONS NCHAR=$nchar;\n\tFORMAT DATATYPE = STANDARD GAP = - INTERLEAVE = yes MISSING = ? SYMBOLS = \"  0 1\";\n\tMATRIX\n";
	}
	$i++;
	print OUTPUT "\n[Matrix generated from $TreeFileName and $SpecTreeFile]\n";
	foreach my $keys (sort keys %Matrix) {
		print OUTPUT "$keys          @{$Matrix{$keys}}\n"; 
	}
	foreach my $kys (sort {$a <=> $b} keys %Comments) {
		$tre = $nchar + $kys;
		print OUTPUT "[$tre: is coded from gene tree $TreeFileName protein(s) $Comments{$kys}\t\t\t\tSistertaxa: $Sisters{$kys}]\n";
	}
	$nchar = $tre;
	close OUTPUT;
	undef %Matrix;
	return ($i,$nchar);
}

sub finish_Matrix {
	my $file = shift;
	my $nchar = shift;
	open MATRIX, "<$file" || die "creation of output file failed: $!\n";
	while (<MATRIX>) {
		if ($_ =~ /(\tDIMENSIONS\sNCHAR.*)/) {
			$hit = $1;
			$hit =~ s/NCHAR=0;/NCHAR=$nchar;\n/;
			push (@outtext, $hit);
		} else {
			push (@outtext, $_);
		}
	}
	close MATRIX;
	open NEWMATRIX, ">$file" || die "creation of output file failed: $!\n";
	print NEWMATRIX "@outtext\nEND;\n";
	close NEWMATRIX;

}
              
sub is_monotypic {
	my $noderef = shift;
	my $node = $$noderef;
	my @gyereklist = ();
	foreach my $gy (my @gyerekek = @{$node->get_terminals}) {
		$gyname = $gy->id;
		if ($gyname =~ /$regex/i){ 
			$dID = $1;
			$dID =~ s/_//; 
			push(@gyereklist, $dID);
		}
	}
	my @uns = keys %{{ map{$_=>1}@gyereklist}};
	my $Arr_length = scalar @uns;
	return(\$Arr_length);
}

sub ancestry {
	my $back_ref = shift;
	my $treeref = shift;
	my @backtrack = @$back_ref;
	my $tree = $$treeref;
	my @ofspecies = @sisters = @uns = ();
	my $node = $backtrack[$#backtrack]; #node is the internal id of the last node in the backtrack array
	my $rootnode = $tree->get_root();
	if ($rootnode ne $SpecNodes_h{$node}) { 
		my $parent = $SpecNodes_h{$node}->get_parent;
#		my @kids = @{$parent->get_children};
		foreach (my @kids = @{$parent->get_children}) {
			next if ($_  eq $SpecNodes_h{$node});
			foreach $kids (my @kids = @{$_->get_terminals}) {
					my $ezt = $kids->id;
					if ($ezt =~ /$regex/i){ 
						$nodeID = $1;
						$nodeID =~ s/_//;
					}
					push(@ofspecies, $nodeID);
				}
		}
	} else { # if rootnode is the backtrack node -> in this case there is no parent
		foreach (my @kids = @{$SpecNodes_h{$node}->get_children}) {
			my $lastbutone = $backtrack[$#backtrack-1];
			next if ($_ ->is_terminal or $_  eq $lastbutone); #this handles the case when the last backtrack node is the root node and thus it has no parent. 
			foreach $kids (my @kids = @{$_->get_terminals}) {
					my $ezt = $kids->id;
					if ($ezt =~ /$regex/i){ 
						$nodeID = $1;
						$nodeID =~ s/_//;
					}
					push(@ofspecies, $nodeID);
				}
		}
	}
	@uns = keys %{{ map{$_=>1}@ofspecies}};
	return(\@uns);
}

sub last_ancestry {
	my $term_ref = shift;
	my $terminal = $$term_ref;
	my $ghj;
	my @ofspecies = @uns = @sisters = ();
	foreach (my @sisters = @{$terminal->get_sisters}) {
		if ($_ ne $terminal) {
			foreach $kids (my @kids = @{$_->get_terminals}) {
				my $ezt = $kids->id;
				if ($ezt =~ /$regex/i){ 
					$nodeID = $1;
					$nodeID =~ s/_//;
				}
				push(@ofspecies, $nodeID);
			}
			@uns = keys %{{ map{$_=>1}@ofspecies}};
		}
	}
	return(\@uns);
}

sub find_root {
	my $in = shift;
	my $intree = $$in;
	my $midpoint = $intree->get_midpoint;
	$intree->reroot($midpoint);
#	foreach my $innode (my @innodes = $intree->get_nodes()) {
#		next if ($innode->is_terminal);
#		my @desclist = @{ $innode->get_children };
#		print "innodes @desclist\n";
#		if ($#desclist == 0) {
#			my $rotnode = $intree->get_root();
#			$intree->set_root_node($innode);
#			my $success = $intree->remove_Node($rotnode);
#			print "suc $success $intree\tt\n"; 
#			my $output = $intree->to_nexus;
#			print "$output";
#			# my $in = new Bio::TreeIO(-format => 'newick', -file => 'tree.tre');
#			$intree->resolve();
#			print "removed the stupid extre node\n";
#		}
#	}
	$intree->resolve;
	return(\$intree);
}

sub mytime {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	$dayOfYear = $daylightSavings;
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	return($theTime);
}

sub small_clusters {
	my $clustfile = shift;
	my $exc = shift;
	my %Excluded = %$exc;
	my $rank = 0;
	open CLUSTER, "<$clustfile";
	while (<CLUSTER>) {
		undef @Matrix{@SpecIDs};	# set an undef value for each species of the species tree in the big matrix
		undef %Comments;
		undef %Sisters;
		$rank++; #this is to count the cluster number
		print " Processing Cluster$rank\n";
		my @Elements = split(' ', $_);
		my $length = scalar @Elements;
		my %ProtPairs = ();
		if ($length < 4) { 
			next if (defined $Excluded{$rank});
			my @proteins = ();
			foreach (@Elements) {
				if ($_ =~ /$regex/i){  #28580|Mellp1
					my $lID = $1;
					$lID =~ s/_//; 
					push(@proteins, $lID);
					$ProtPairs{$lID} = $_;
				}
			}
			@uns = keys %{{ map{$_=>1}@proteins}};
			&code_it(\%ProtPairs,\@uns);
			$unss = scalar @uns;
			if ($unss <= 2) {
				foreach (keys %ProtPairs){
				}
				foreach my $elem (@Elements) {
					my %LastPairs = ();
					my $rID;
					if (!grep {$_ eq $elem} values %ProtPairs) {#!exists $ProtPairs{$_}) {
				 		if ($elem =~ /$regex/i){  
							$rID = $1;
							$rID =~ s/_//; 
						}
				 		$LastPairs{$rID} = $elem;
				 		my @sis = keys %LastPairs;
				 		if ((scalar keys %LastPairs) > 0) {
				 			&code_it(\%LastPairs,\@sis);
				 		}
				 	}
				}
			}
			my $ClustName = 'Cluster' . "$rank";
			($i, $nchar) = &Dumper($ClustName, $OutFileName, $i);

		}
	}
	close CLUSTER;
	return($nchar);
}