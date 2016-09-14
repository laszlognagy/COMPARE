#!/usr/bin/perl

use strict;
use warnings;
use Parallel::ForkManager;
use Getopt::Long;

my ($AliDir, $TreeDir, $Outdir, $SpecTree, $NProc, $restart, $treefile);
my (@Files_done);

#USAGE: ./treefix_wrap.pl --alidir [directory with alignment] --treedir [gene tree dir] --outdir [dir for output tree files] --spectree tree_file --nproc integer [--restart]
# alidir: directory containing alignments
# treedir: dir containing rooted, binary gene trees
# spectree: species tree used for reconciliation
# nproc: number of CPUs
# May need to modify lines 64 & 73 depending on the file extension and naming convention used. Also check the regexp-s in the create_smap sub and modify if needed.
# note that | signs suck at treefix, so they have to be replaced in the input files - see reformat.pl for this purpose

###############Parsing command line arguments############
GetOptions( 'alidir=s' => \$AliDir, #currently unused
			'spectree=s' => \$SpecTree, # 
			'treedir=s' => \$TreeDir, # 
			'nproc:s' => \$NProc, #path and name of bayestraits version
			'restart' => \$restart,
			'outdir=s' => \$Outdir); #currently unused

# Check command line parameters
if (not $SpecTree or not (-e $SpecTree)) { die "ERROR: missing input organismal tree file\n"; }
#if (not $AliDir or not (-d $AliDir)) { die "ERROR: cant find input alignment directory\n"; }
if (not $TreeDir or not (-d $TreeDir)) { die "ERROR: cant find input tree directory\n"; }
if (not $Outdir) { 
	print "You have not specified a directory for results, I will use /tree_files_rearranged\n"; 
	$Outdir = 'tree_files_rearranged';
}
#unless (-d "$Outdir") {
#	mkdir $Outdir or die "Unable to create $Outdir: $!\n";
#}
if (not $NProc) { 
	print "You have not specified the number of parallel processes to run (nproc argument) \n"; 
	$NProc = 1;
}

if ($restart) {
	open RELOAD, "<inferences_completed.tmp";
	while (<RELOAD>) {
		chomp;
		push(@Files_done, $_);
		print "Restarting calculations from the last file $_\n";
	}
}


sleep(5);


#my $SpecTree = "71GA.tre.tree";
#my $NProc = 1;
#my $startdir = ".";
#my $treefile = ();
#chdir "/Files";
my $Phy = new Parallel::ForkManager($NProc);
opendir(TREEDIR, $TreeDir) || die "can't opendir alignment dir: $!";
while($treefile = readdir(TREEDIR)) {
	if ($treefile =~ /.trimmed.tre$/) { #if there are failed tree inferences originally in raxml as in my first case, it will not cause problems here, even if the alignments are in the directory
		if ($restart) { # Skip trees that have alrady been reconciled
			next if (grep{$_ eq $treefile}@Files_done);
		}
		print "file being reconciled is $treefile\n";
	    my $pid = $Phy->start and next;   # Forks and returns the pid for the child:
			&create_smap($TreeDir, $treefile);
			my $logfile = $treefile . '.treefix.log';
#			&parse_xml($alifile); (ez az excludefile krealasahoz kellene?)
			my $task = 'treefix' . ' -s ' . $SpecTree . " -S smapfile$treefile.smap" . ' -A .out.best.fas_trimmed -o .out.best.fas_trimmed.tre -n .out.treefix.tree -V1 -e \'-m PROTGAMMAWAG\''  . ' -l ' . "$logfile $treefile";
			my $return = 0;
			#system ($task);
			print "$task\n";
	    if ($return == 0) {
	    	open READYTREES, ">>inferences_completed.tmp" || die "Cant create temporary file: $!";
			print READYTREES "$treefile\n";
			close READYTREES;
		} else {
	    	open FAILTREES, ">>inferences_failed.tmp" || die "Cant create temporary file: $!";
			print FAILTREES "$treefile\n";
			close FAILTREES;			
		}
	 #   unlink ("smapfile.smap");
	    $Phy->finish; # Terminates the child process
	}
}
$Phy->wait_all_children;
print "A list of successfully reconciled trees have been printed to inferences_completed.tmp\n";
#unlink("*.smap");

exit;

sub create_smap {
	my $dir = shift;
	my $file = shift;
	$file = './' . $dir . '/' . $file;
	open TREE, "<$file" || die "ezitten $!";
	my @matches = ();
	my $str = ();
		print "treefile is $file\n";
	while (<TREE>) {
		if ($_ =~ /(\d*_\w*\d?):/) {
#			print "match is $1\n";
		}
		push (@matches,$&) while($_ =~ /(\d+_\w+\d?)/g);
#		print "matches @matches\n";
#		@matches = $str =~ /,|\((\d*_\w*):/g;
			open SMAP, ">>smapfile$treefile.smap" || die "opening file for smap failed: $!\n";
			foreach (@matches) {
				s/://;
				my $seqn = $_;
				s/\d*\_//;
				print SMAP "$seqn\t$_\n";
			}
			close SMAP;
	}
	close TREE;
}