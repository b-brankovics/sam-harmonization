#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A wrapper for NCBI's BLAST that produces valid SAM files as output.\n" .
    "\tOnly blastn can be run using this wrapper, the others are difficult\n" .
    "\tto translate to SAM. If <db> has not be indexed by makeblastdb, then\n" .
    "\tthat will be done during the run.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] <db> <query>\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-a | --all\n\t\tKeep all query sequences if they have no hits.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

my @arg;
my $keep_all;
for (@ARGV) {
    if (/^--?a(ll)?$/) {
	$keep_all++;
    } else {
	push @arg, $_;
    }
}

@ARGV = @arg;

# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my ($ref, $query) = @ARGV;

my $blastdb = $ENV{'BLASTDB'};

#print STDERR "BLASTDB=$blastdb\n";

unless ($ref && $query) {
    biointbasics::print_help(\@ARGV, $info, "ERROR: specify at least two input FASTA files.\n");
} else {
    if (! -e $ref || -z $ref) {
	if ( -e $ref . '.nal' || -e $ref . '.nin') {
	    # It is a DB in the current folder
	} elsif ($blastdb && (-e $blastdb . $ref . '.nal' || -e $blastdb . $ref . '.nin') ) {
	    print STDERR "INFO: There is no '$ref' file in the current folder, but there is BLAST DB in '$blastdb'\n";
	} else {
	    biointbasics::print_help(\@ARGV, $info, "ERROR: '$ref' does not exist or it is empty.\n");
	}
    }
    if (! -e $query || -z $query) {
	biointbasics::print_help(\@ARGV, $info, "ERROR: '$query' does not exist or it is empty.\n");
    }
}

my $makeblastdb = 'makeblastdb';
my $versiontext = `$makeblastdb -version`;
$versiontext =~ /^\S+: (\S+)/;
my $version = $1;
my $noindex;
for (qw/.nhr .nin .nog .nsd .nsi .nsq/) {
    my $file = $ref . $_;
    if (! -e $file || -z $file) {
	$noindex++;
    }
}
if ($noindex && -e $ref . '.nal' && ! -z $ref	. '.nal') {
    # Skip indexing if there is an alias file ('.nal')
    $noindex = undef;
}
if ($noindex && -e $ref . '.nin' && ! -z $ref	. '.nin') {
    # Skip indexing if there is an index file ('.nin'); update_blastdb.pl does not create the same extensions like when you create a DB from FASTA
    $noindex = undef;
}
if ($noindex && $blastdb && (-e $blastdb . $ref . '.nal' || -e $blastdb . $ref . '.nin' ) ) {
    $noindex = undef;
}
if ($noindex) {
    my $mkdb_cmd = $makeblastdb . " -dbtype nucl -parse_seqids -in $ref";
    `$mkdb_cmd`;
}
my $blastn = 'blastn';
my $outfmt = 6;
$outfmt = 7 if $keep_all;
my $cmd = $blastn . " -db $ref -parse_deflines -query  $query -outfmt " . '"' . $outfmt . ' sseqid slen sstart send qseqid qlen qstart qend sstrand mismatch gaps btop sseq qseq score bitscore evalue"';
#print STDERR "$cmd (BLASTDB=$blastdb)\n";

my $output = `$cmd`;
#$output =~ s/\n([A-Za-z]+)/$1/g;
#open my $input, "$cmd |";
biointsam::parse_blast($output, $cmd, $version);



sub reverse_seq {
    # Reverse complements the sequences
    my ($seq) = @_;
    # Reverse the sequnce
    my $complement = reverse $seq;
    # Complement the sequence
    $complement =~ tr/ACGTacgtWwMmRrSsKkYyBbVvDdHh/TGCAtgcaWwKkYySsMmRrVvBbHhDd/;
    return $complement;
}
