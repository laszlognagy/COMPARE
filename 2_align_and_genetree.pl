#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;

my (@Clusters_done, @Clusters_aligned, @Trees_done, @Unverified_clusters, $NProc, $inputdir, $AliDir, $exclude, $TreeDir, $Outgroups, $Aligner, $Treeference, $restart, $alifile, $Model);

#USAGE: ./runner.pl --input_dir <indir> [--ali_dir <alignmentdir>] [--tree_dir <treedir>] --ncpu <integer> [--aligner <program_name>] [--ml_prog <program name and path>] [--exclude] [--restart]
# input_dir: directory containing unaligned protein fasta files
# ali_dir: directory where completed alignments are written (directory must exist)
# tree_dir: directory where completed gene trees are written (directory must exist)
# ncpu: number of CPU-s to use
# ml_prog: path to the tree inference program of choice
# exclude: if specified clusters listed in the file on line 56 (Clusters_with_no_significant_similarity.txt) will be skipped. Useful to e.g. exclude clusters that contain protein that show no similarity to anything else.
# restart: restart if needed. If requested, will read the list of completed tasks from ./TempFiles/Clusts_aligned.tmp. Make sure permissions are set to allow creating this dir.
# May need to modify lines 84, 95, 96, 111, 121, 134 depending on the type and path to the MSA and gene tree inference software.
#WARNING: If using prank for alignment, it is advisable to check manually if no alignments failed in the range of the 50 biggest clusters. PRANK may die without error message in some cases, mostly with clusters >200 proteins. In these cases, the missing alignments have to be added manually. An automatic checking will be included in the future versions.

#parsing command line arguments
GetOptions( 'input_dir:s' => \$inputdir, #Working Directory
			'ali_dir:s' => \$AliDir, #alignment output Directory
			'tree_dir:s' => \$TreeDir, #tree output Directory
			'ncpu=i' => \$NProc, #Number of processes to run in parallel
#			'outgr=s' => \$Outgroups, #text file containing outgroup information
			'aligner:s' => \$Aligner, #Aligner is the path and name of the alignment program to use - default will be prank in a folder in the package. Can be substituted with anything, especially if speed is of concern. Its advisable to install exonerate, however when the default prank version is chosen. Also, if you want to use the +F option in PRANK, you will need to modify the script yourself.
			'ml_prog:s' => \$Treeference, #Treeference if the path and name of the raxml implementation  - default will be prank in a folder in the package
			'exclude' => \$exclude, #whether to exclude unverified clusters or not
			'restart' => \$restart); # whether to restart or not
#print "first $ClustFile $FullNames $MinSize $restart $outfilder $phylodir $allgenomes\n";

# Check command line parameters
#if (not $Outgroups or not (-e $Outgroups)){ die "ERROR: missing input outgroup file\n"; }
if (not $inputdir){ 
	$inputdir = "Cluster_files";
	print "WARNING: You have not specified a working directory, so I assume your input fasta files are in the current folder\n";}
if (not $AliDir){ 
	print "You have not specified an alignment directory, a folder called Alignment_files will be created for you\n"; 
	$AliDir = 'Alignment_files';}
unless (-e $AliDir) {
	mkdir $AliDir;}
if (not $TreeDir){ 
	print "You have not specified a tree output directory, a folder called Tree_files will be created for you\n"; 
	$TreeDir = 'Tree_files';}
unless(-e $TreeDir) {
	mkdir $TreeDir;}
if (not $Aligner){ 
	print "You have not specified a path to the aligner, PRANK in the ./executables folder will be used\n"; 
	$Aligner = '/executables/prank';}
if (not $Treeference){ 
	print "You have not specified a path to the tree inference program, RAxML in the ./executables folder will be used\n"; 
	$Treeference = '/executables/raxml';}

if ($exclude) {
	open EXCFILE, "<./TempFiles/Clusters_with_no_significant_similarity.txt" || die "Could not open ";
	while (<EXCFILE>) {
		chomp;
		push(@Unverified_clusters, $_);
	}
	close EXCFILE;
}

if ($restart) {
	open RELOAD, "<./TempFiles/Clusts_aligned.tmp" or die "Cannot restart because a list of already completed clusters is missing";
	while (<RELOAD>) {
		chomp;
		push(@Clusters_done, $_);
	}
	close RELOAD;
	if (-e "./TempFiles/Tree_inferences_completed.tmp") {
		open RELOOD, "<./TempFiles/Tree_inferences_completed.tmp";
		while (<RELOOD>) {
			chomp;
			push(@Trees_done, $_);
		}
	close RELOOD;
	}
}

my $Ali = new Parallel::ForkManager($NProc);
opendir(INDIR, $inputdir) || die "can't opendir $inputdir: $!";
while(my $file = readdir(INDIR)) {
	if ($file =~ /^Cluster.*/) {
		if ($restart) { # Skip clusters that have alrady been aligned
			if (grep{$_ eq $file}@Clusters_done) {
				print "$file already aligned, skipping...\n";
				next;
			}
		}
		if ($exclude) { #Skip clusters which are unverified
			next if (grep{$_ eq $file}@Unverified_clusters);
		}
	    my $pid = $Ali->start and next;   # Forks and returns the pid for the child:
			my $outfile = $file . 'aligned.out';
			my $indeptask = "/opt/prank-msa/bin/prank" . ' -d=' . './' . $inputdir . '/' . $file . ' -o=' . './' . $AliDir . '/' . $outfile . ' -quiet -showxml -support'; # -showxml -support to get posterior scores
			# lehet hogy az exonerate beszarik multiple processes alatt????? csekkolni megegyszer es betenni a -noanchors argumentet ha igy van
			print "Alignment of $file started\n";
			print "task $indeptask\n";

#			system ($indeptask); # valahogy prank return value-kat megnezni es ha failed az alignment akkor extra kommenttel kiirni a clucts aligned fajlba
		open READYS, ">>./TempFiles/Clusts_aligned.tmp" || die "Cant create temporary file: $!";
		print READYS "$file\n";
		close READYS;
	    $Ali->finish; # Terminates the child process
	}
}

opendir ALDIR, $AliDir;
while(my $file = readdir(ALDIR)) {
	if ($file =~ /.*out.1.fas$|.*out.1.xml$/) { #CHECK THIS
		unlink $file;
	}
}

#valahogy ckeckolni hogy keszen van-e mindegyik ALIGNMENT???

my $Phy = new Parallel::ForkManager($NProc);
opendir(ALIDIR, $AliDir) || die "can't opendir alignment dir: $!";
while(my $alifile = readdir(ALIDIR)) {
	if ($alifile =~ /fastaaligned.out.2.fas$/) {
		if ($restart) { # Skip clusters for which trees have alrady been inferred
			if (grep{$_ eq $alifile}@Trees_done) {
				print "ML tree for $alifile already exists....\n";
				next;
			}
		}
		print "file is $alifile\n";
		my ($filename, $size) = &fasta_to_phylip($alifile, $AliDir);
		if ($size < 50) {$Model = 'PROTGAMMAWAG'} else {$Model = 'PROTCATWAG';}
	    my $pid = $Phy->start and next;   # Forks and returns the pid for the child:
			my $outfile = $alifile . '.tree.out.tre';
#			&parse_xml($alifile);
			my $indeptask = "/opt/RAxML-7.2.8-ALPHA/bin/raxmlHPC" . ' -s ' . $filename . ' -n ' . $outfile . ' -m ' . $Model; # . ' -E excludefile'
			print "Maximum Likelihood estimation for of $filename started\n$indeptask\n";
#			system ($indeptask);
	    open READYTREES, ">>./TempFiles/Tree_inferences_completed.tmp" || die "Cant create temporary file: $!";
		print READYTREES "$alifile\n";
		close READYTREES;
	    rename "RAxML_result.$outfile", "./$TreeDir/$outfile";
	    unlink ("RAxML_bestTree.$outfile", "RAxML_info.$outfile", "RAxML_log.$outfile", "RAxML_parsimonyTree.$outfile");
	    $Phy->finish; # Terminates the child process
	}
}
$Phy->wait_all_children;
print "A list of successfully aligned clusters and inferred trees can be found in the files Clusts_aligned.tmp and Tree_inferences_completed.tmp\n";
exit;

###########Subs###########
sub fasta_to_phylip {
	  # hash lookup of sequences based on taxon names. for each input fasta file a separate hash is created
	my $onematrix = shift;
	my $alidir = shift;
	my %Matrix = ();
	my $path = './' . $AliDir . '/' . $onematrix;
	print "path $path\n";
	local $/ = '>';
	open (FASTA, "<$path") || die "error2: $!";
	while (<FASTA>) {
		s/^>//g; # strip out '>' from beginning
		s/>$//g; # and end of line
		next if !length($_);           # ignore empty lines
		my ($header_info) = /^(.*)\n/; # capture the header
		s/^(.*)\n//;                   # and strip it out
#    	my @rec = split /\|/, $header_info;
		s/\n//mg;                      # join the sequence strings
#   	push @rec, $_;                 # make sequence the last element
#    	push @seq, \@rec;              # push into array of sequences
		$Matrix{$header_info} = $_;         # or store in a hash lookup
	}
	close FASTA;
#	chdir $alidir;
	my $phyfile = './' . $alidir . '/' . $onematrix . '.phy';
	my $charnum = split(//, $Matrix{(keys %Matrix)[rand keys %Matrix]});
	my @taxlist = keys %Matrix;
	my $taxnum = scalar @taxlist;
	open PHYLIP, ">$phyfile" || die "Could not create temporary phylip file in $alidir folder, please check permissions\n";
	print PHYLIP "$taxnum $charnum\n" or die "cannot write to phylip output file $phyfile\n";
	foreach (keys %Matrix) {
		print PHYLIP "$_         $Matrix{$_}\n";
	}
	close PHYLIP;
	return($phyfile, $taxnum);
}

#sub parse_xml {
#	my $infile = shift;
#	$infile = s/.out.2.fas/.out.2.xml/;
#	my ($counter, @probrec, %xml);
#	if ($exclude) {
#		open XML, ">$infile";
#		local $/ = '  </probability>';
#		while (<XML>) {
#			s/^  <\/probability>//g; # strip out '>' from beginning
#			s/  <\/probability>$//g; # and end of line
#			if ($_ =~ /<probability id="2">/) {
#				s/  <probability id="2">\n    //;
#				$counter++;
#				$probrec = length;
#				$xml{$counter} = $_;
#			}	
#		} #here the hash should be ready
#		loop through individual sites
#		for (0 .. $probrec) { #something like this
#			get the sum or mean of the postprobs for all nodes for that site
#		}
#	} else {
#		open EXCL, ">excludefile";
#		print EXCL "1-1";
#		close EXCL;}
#
#
#
#}