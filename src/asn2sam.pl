#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;

#===DESCRIPTION=================================================================
# A tool to remove duplicate sequence from FASTA format files and
#  print the groups to STDERR

my $description = 
    "Description:\n\t" .
    "A tool to convert BLAST archive (ASN.1) files to SAM format and\n" .
    "\tprint to STDOUT\n";
my $usage = 
    "Usage:\n\t$0 [-h | --help] [-a | --all]  <ASN.1 file>\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-a | --all\n\t\tPrint sequences as unmapped that had no hits in the BLAST\n" .
    "\t\t(By default these sequences are not reported.)\n" .
    "\n";

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

# Print help if needed
&print_help(\@arg);

if (@arg != 1 ) {
    my $msg = "ERROR: There has to be exactly one ASN.1 file given as input\n";
    &print_help(['-h'], $msg);
} elsif (! -e $arg[0]) {
    my $msg = "ERROR: '$arg[0]' does not exist!\n";
    &print_help(['-h'], $msg);
}    
my ($archive) = @arg;

my $blastfmt = 'blast_formatter';
my $versiontext = `$blastfmt -version`;
$versiontext =~ /^\S+: (\S+)/;
my $version = $1;
my $outfmt = 6;
$outfmt = 7 if $keep_all;
my $cmd = $blastfmt . " -archive $archive -parse_deflines -outfmt " . '"' . $outfmt . ' sseqid slen sstart send qseqid qlen qstart qend sstrand mismatch gaps btop sseq qseq score bitscore evalue"';

my $output = `$cmd`;
unless (${^CHILD_ERROR_NATIVE} == 0) {
    die "ERROR: '$blastfmt' has crashed. Check if you are using the correct input and the BLAST DB used to generate the ASN is where it has to be.\n" 
}
biointsam::parse_blast($output, $cmd, $version);

#===SUBROUTINES=================================================================

sub print_help {
    # Print out the usage to STDERR
    # Takes in the @ARGV as input
    my ($args, $error) = @_;
    $error = "" unless $error;
    for (@$args) {
	if (/^-?-h(elp)?$/) {
	    die $error . "$usage\n$description\n$options";
	}
    }
}
