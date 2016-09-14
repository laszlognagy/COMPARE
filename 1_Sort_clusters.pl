#!/usr/bin/perl

# USAGE: sort_cluster.pl --cf <cluster_file> --fn <filename> --wd <working directiry> --merged <filename> --clud <directory for clusters> --min <integer> --phyd <directory for singlecopy clusters> --restart
# TAKES AS INPUT: 
#	- a file specifying the full names of each species (option --fn) Only needed if you let the script run blast searches.; 
#	- a cluster file in which prot ID-s in one line corresponds to one cluster (option --cf);
#	- optionally a working directory if different from the current directory (option --wd). the working directory should contain the input fasta files. input fasta files should end with "input.fasta";
#	- optionally a directory to output assembled clusters (option --clud), default setting is Cluster_files;
#	- optionally a filename to save the merged genomes (option --merged), default is allgenomes.fasta. Note this feature is pointless if the genomes are already merged in one file
#	- optionally a directory where singlecopy clusters which are potentially useful in inferring species phylogeny can be sorted (option --phyd);
#	- optionally an integer value used as the maximum cluster size not to subject to BLAST-based verification (option --min). Default value is 5. If you want to avoid blast set --min to 0;
#	- if the --restart option is specified the script continues its job from where it has been stopped
# PRODUCES: 
#	- one fasta file for each cluster in the directory specified on the command line; 
#	- a merged fasta file containing proteins from all input genomes;
#	- a file listing clusters for which no significant similarity could be found by BLAST (if blast is allowed);
#	- on demand: a directory containing clusters suitable for species tree inference;

# !!!!NOTE the regexp use throughout the script is \w*\_(\w*), which captures protein names in the form 32439845_Lacbi1, i.e., protein ID + underscope + speciesname. Modify regexp if needed.

#use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Tools::Run::RemoteBlast;
use Getopt::Long;

my %Matrix;
my ($FullNames, $workdir, $outfilder, $ClustFile, $MinSize, $allgenomes, $phylodir, $restart);

#parsing command line arguments
GetOptions( 'fn=s' => \$FullNames, #translation table
			'wd:s' => \$workdir, #Working Directory
			'clud:s' => \$outfilder, #ClusterDirectory
			'cf=s' => \$ClustFile, #ClusterFile
			'min:i' => \$MinSize, #minimum size of cluster to subject to blast-based check
			'merged:s' => \$allgenomes, #merged genomes fasta file
			'phyd:s' => \$phylodir, #single copy directoiry
			'restart' => \$restart); # whether to restart or not
#print "first $ClustFile $FullNames $MinSize $restart $outfilder $phylodir $allgenomes\n";

# Check command line parameters
if (not $ClustFile or not (-e $ClustFile)){ die "ERROR: missing input cluster file\n"; }
if (not $FullNames or not (-e $FullNames)){ die "ERROR: missing input species translation table\n"; }
if (not $workdir){ 
	$workdir = ".";
	print "WARNING: You have not specified a working directory, so I assume your input fasta files are in the current folder\n";
}
if (not $outfilder){ 
	print "You have not specified an output directory, a folder called Cluster_files will be created for you\n"; 
	$outfilder = 'Cluster_files';
}
if (not $allgenomes){ 
	print "You have not specified a name for the merged fasta file, merged genomes will be saved to Allgenomes.fasta\n"; 
	$allgenomes = 'allgenomes.fasta';
}
if (!$MinSize) { $MinSize = 3; }

#print "workdir $workdir\n";
#print "$ClustFile $FullNames $MinSize $outfilder $phylodir $allgenomes\n";

#my $Clustercounter = 0;
#$restart = 'no'; # this should be loaded as a command line parameter
#my @Cluster = ();
#my $ClustFile = '/Users/laszlonagy/Documents/pipe/cluster/SAP_31_clustering2/MCL_2/mcl.cls_2_orig'; #Or read from command line option
#my $workdir = '/Users/laszlonagy/Documents/pipe/phylo/raw_proteins'; #or simply say "."
#my $outfilder = "clusterfiles";
#my $phylodir = "singlecopy_clusters";


if ($restart) { # the checkpointing option should be put into the other oroginal version too
	opendir(CHECKPOINT, $outfilder);
	while(my $r = readdir(CHECKPOINT)) {
		if ($r =~ /^Cluster.*/) {
			push(@Clusters_ready, $r); # clusters ready contains the filenames of already done clusters
		}
	}
} else {
	mkdir $outfilder || die "cannot create directory: $!";
	}

if (defined $phylodir) { # phylodir should contain the name of the directory in which the singlecopy genes would be dumped, ha nem kerjuk hogy ezt csinalja akkor hylodir undef kell hogy legyen
	mkdir $phylodir  || die "cannot create directory: $!";
}

#Merge Fasta Files
opendir(WORKING, $workdir) || die "can't opendir $workdir: $!";
while(my $file = readdir(WORKING)) { # this part will make the script work with perl 5.11 or newer
	if ($file =~ /.*input.fasta/) { #all input genome fasta files should end with "input.fasta"
		print "reading $file.....\n";
		$path = "$workdir" . '/' . "$file";
		$Matref = &fastaparser($path);
		%Mat = %$Matref;
	}
} #at the end of while loop, we will have a %Mat of all the proteins from all input genomes

#printing merged fasta file
print "Just in case you needed, proteins from all input genomes merged together has been written to $allgenomes\n";
open MERGEDFASTA, ">$allgenomes" || die "error: $!";
foreach (keys %Mat) {
	print MERGEDFASTA ">$_\n$Mat{$_}\n";
}
close MERGEDFASTA;

#Get Full taxon names from provided input file
my %FullNames = ();
open NAMESFILE, "$FullNames";
while (<NAMESFILE>) {
	chomp;
	my @Record = split(/\t/, $_);
	$FullNames{$Record[1]} = $Record[0]; # shange this and then change input file format???
}

#Creating TempFiles folder for all useful and temporary files during the analysis
mkdir "TempFiles" || die "cannot create directory: $!";


#Write clusters to separate files and check clusters by BLAST
open CLUSTFILE, "$ClustFile"; #assumes each line of the cluster file contains members of one cluster separated by tabs
while (<CLUSTFILE>) {
	chomp;
	my @Species = @Clustsize = @SpecNames = ();
	my @Cluster_entries = split("\t", $_);
#	$Clustarr_ref = &Search(\@Cluster_entries);
	my $Clusterfile_out = 'Cluster' . $. . '.fasta';
	if ($restart && grep {$_ eq $Clusterfile_out} @Clusters_ready) {
		print "Cluster $Clusterfile_out already done, skipping.....\n";
		next;
	}
	print "cluster $. assembled....\n";
	&Dumper(\@Cluster_entries, \$Clusterfile_out, \$outfilder); # print to file all clusters but produce list of those with no significant match
	if ($phylodir) { # do this part only when phylodir is defined, i.e. we asked for sorting of phylo-useful clusters
		foreach (@Cluster_entries) {
				$_ =~ /\w*\_(\w*)/; # grab species name jgi|Copci1|1|CC1G_00004T0
#				print "grabbed for comparison is $1\n";
				push(@SpecNames, $1);
			}
		@Clustsize = keys %{{ map{$_=>1}@SpecNames}};
		if ((scalar @Clustsize) == (scalar @Cluster_entries)) {
			print "Found cluster of single-copy genes, printing to directory $phylodir\n";
			&Dumper(\@Cluster_entries, \$Clusterfile_out, \$phylodir);
		}
	}
#	my @Cluster = @$Clustarr_ref;
	if (@Cluster_entries <= $MinSize) {
		my $real_hit_count = 0;
		print "Found cluster smaller than specified minimum size\n";
		foreach (@Cluster_entries) {
			$_ =~ /\w*\_(\w*)/; # grab species name jgi|Copci1|1|CC1G_00004T0
#			print "grabbed is $1\n";
			push(@Species, $1);
		}
		my @spec_labs_uniq=keys %{{ map{$_=>1}@Species}};
#		print "uniqs is @spec_labs_uniq\n";
		if ((scalar @spec_labs_uniq) > 1) { #Means non-monotypic cluster
			print "Cluster $Clusterfile_out is not monotypic, skipping...\n";
			next;
		}
		if ((scalar @spec_labs_uniq) == 1) { #if cluster monotypic
#			print "innen blast\n"; #osszerakni a
			$FileNamm = './' . $outfilder . '/' . $Clusterfile_out;
			my $stream = Bio::SeqIO->new(-file=> $FileNamm,-format=>'fasta');
#			print "read succesd\n";
			while (my $seq = $stream->next_seq) {
				my $remote_blast = Bio::Tools::Run::RemoteBlast->new (
      												     -prog => 'blastp',-data => 'nr',-expect => '1e-10' ); #, -readmethod => 'xml'
				my $rc = 
#				print "remoteblst: $remote_blast\n";
				my $r = $remote_blast->submit_blast($seq); #first ensure oinly single seq files are submitted
				print "Waiting for BLAST results.....\n";
				while ( my @rids = $remote_blast->each_rid ) { #analyzing blast results
#					print "miaza rid: @rids\n";
 					foreach my $rid ( @rids ) {
 #						print "remote blast id is $rid\n";
 						my $rc = $remote_blast->retrieve_blast($rid); #rc contains either a constant or Bio::SearchIO object
#						print "erce is $rc\n";
						if( !ref($rc) ) {
 							if( $rc < 0 ) {
							$remote_blast->remove_rid($rid);
							}
 							sleep 5;
						} else {
							my $result = $rc->next_result();
#							print "result is $result\n";
							if (!$result) { #ha result ures -> no matches found
								print "it seems no significant hit for this query was found\n";
								$real_hit_count++;
							} else {	###IMPORTANT this else should be placed in the original version too
							$remote_blast->remove_rid($rid);
							my $quer_n = $result->query_name();
#							print "\nQuery Name: ", $quer_n, "\n";
							if ($quer_n =~ /\w*\_(\w*)/) {$quer_specn = $1;}
#							print "query specname $quer_specn \n";
							while ( my $hit = $result->next_hit ) {
								my $name = $hit->name;
								my $descr = $hit->description;
	#							print "\thit name is $name\t$descr\n";
								if ($descr =~ /$FullNames{$quer_specn}/) {
#									print "Hit species $FullNames{$quer_specn} and query species\t\t$descr same\n";
									$real_hit_count++;
								} else {
									$real_hit_count = 0;
#									print "valszeg nom atch\n";
									last;
								}
							}
							}
						}
					}
				}
#				print "osszeg is $real_hit_count\n";
			}	
	
		}
		if ($real_hit_count > 0) { #probably save this into a separate working folder?
			print "Found bad cluster $Clusterfile_out, printing it to file of unreliable clusters: ./TempFiles/Clusters_with_no_significant_similarity.txt\n";
			open BADCLUSTERS, ">>./TempFiles/Clusters_with_no_significant_similarity.txt" || die "Could not create output file: $!";
			print BADCLUSTERS "$Clusterfile_out\n";
			close BADCLUSTERS;
		}
	} #if<5
}

##############Subs#############

sub fastaparser {
	  # hash lookup of sequences based on taxon names. for each input fasta file a separate hash is created
	my $onematrix = shift;
	local $/ = '>';
	open (FASTA, $onematrix) || die "error2: $!";
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
	return (\%Matrix);
}

sub Search {
	my $s = shift;
	my ($SpecName, $ProtID);
	my @Array_of_entries_in_cluster = ();
	my @Cluster = @$s; # 53570|Trave1
	foreach $element (@Cluster) {
#		my @elem = split(/\|/, $element); #ASSUMES format!!!
#		my $query = "jgi|" . $elem[1] . "|" . $elem[0]; #query is the ProtID|Specname formatted entry coming from the clustering file
		$element = quotemeta($element);
		foreach $hat (keys %Mat) {
			$hat = quotemeta($hat);
#			if ($hat =~ /jgi\|(\w*)\|(\w*)\|?.*/) {
#				$SpecName = $1;
#				$ProtID = $2;
#			}
#			$r = 'jgi\|' . $SpecName . '\|' . $ProtID;
			if ($element eq $hat) {
#				print "$query matched a key $hat in mat\n";
				push(@Array_of_entries_in_cluster, $hat);
			}
		}
	}
	return(\@Array_of_entries_in_cluster);
}

# if (exists %Mat{/$query/})

sub Dumper {
	my $ref = shift;
	my $cref = shift;
	my $filderref = shift;
	my @Cluster_and_seqs = @$ref;
	my $Clusterfile_out = $$cref;
	my $outfilder = $$filderref;
	chdir $outfilder || die "cannot cd into directory: $!";
	open CLUSTER, ">$Clusterfile_out" || die "error: $!";
	foreach (@Cluster_and_seqs) {
		print CLUSTER ">$_\n$Mat{$_}\n";
	}
	close CLUSTER;
	print "Cluster has been written to $Clusterfile_out\n";
	chdir "..";
#	return ($Clusterfile_out);
}

