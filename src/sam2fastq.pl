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
    "A tool convert to SAM format files to FASTQ format printed on STDOUT and STDERR.\n" . 
    "\tDisclaimer: secondary alignments are skipped, and proper pair order is not validated.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] [SAM file] >reads_R1.fq 2>reads_R2.fq\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my %ref;
my %query;

my $minlen;
my $minaln;
my $minsim;

while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    unless (%hit) {
	#print "$_\n";
	next;
    }
    if ($hit{'FLAG'} & 16) {
	$hit{'SEQ'} =  biointsam::reverse_seq( $hit{'SEQ'} );
	$hit{'QUAL'} = reverse $hit{'QUAL'};
    }
    if ($hit{'QUAL'} eq '*') {
	$hit{'QUAL'} = 'I' x length( $hit{'SEQ'} );
    }
    # Skip secondary alignements
    next if 256 & $hit{'FLAG'};
    # Split forward and reverse between STDOUT and STDERR
    if (64 & $hit{'FLAG'}) {
	print STDERR '@' . $hit{'QNAME'} . "\n", $hit{'SEQ'} . "\n", "+\n", $hit{'QUAL'} . "\n";
    } else {
	print '@' . $hit{'QNAME'} . "\n", $hit{'SEQ'} . "\n", "+\n", $hit{'QUAL'} . "\n";
    }
}
