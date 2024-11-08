#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

use Bio::SeqIO;
use Bio::Seq;

my $programname = "sam-update-cigar.pl";
my $version = "1.1";
my $cmd = join(" ", $programname, @ARGV);

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to update the CIGAR portion of a SAM file to the basic version of CIGAR.\n";
my $usage = 
    "Usage:\n\t$0 [-h] <SAM file>\n";
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

my ($sam) = @ARGV;

# Print help if needed
biointbasics::print_help(\@ARGV, $info);


my %ref;
my %query;
my $hits;
my $header = "true";
my $previousprogram = "";
while(<>) {
    #print "$_";
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit);
    # Header lines contain no hits
    unless (%hit) {
	# Print header to maintain a valid SAM output
	print "$_\n";
	if (/^\@PG\tID:(\S+)/) {
	    $previousprogram = $1;
	}
	next;
    }
    if ($header) { # First line after the header section
	my $text = "\@PG\tID:$programname\tPN:$programname";
	$text .= "\tPP:$previousprogram" if $previousprogram;
	print $text . "\tVN:$version\tCL:$cmd\n";
	$header = undef;
    }

    # Nothing to do if there is no hit for the query sequence, so skip it
    unless($hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*") {

	# Update CIGAR if needed
	my $len = 0;
	my $type;
	my $cigar = "";
	my $original = $hit{'CIGAR'};
	while (length $original) {
	    $original =~ s/^(\d+)//;
	    my $l = $1;
	    $original =~ s/^(\D)//;
	    my $t = $1;
	    if ($t eq '=' || $t eq 'X') {
		$t = 'M';
	    }
	    if ($type) {
		if ($t eq $type) {
		    $len += $l;
		} else {
		    $cigar .= $len . $type;
		    $len = $l;
		    $type = $t;
		}
	    } else {
		$len = $l;
		$type = $t;
	    }
	}
	$cigar .= $len . $type;
	$hit{'CIGAR'} = $cigar;
    }

    print biointsam::sam_string(\%hit), "\n";
}


