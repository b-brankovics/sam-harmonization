#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

my $w = 100;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to display refernce coverage for hits in a SAM file for exploratory purposes.\n";
my $usage = 
    "Usage:\n\t$0 [-h] <SAM file>\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-w=<int> | --width=<int>\n\t\tSpecify the width of the coverage plots to be displayed (default is $w).\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================


# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my @keep;
for (@ARGV) {
    if (/^--?w(?:idth)?=(\d+)$/) {
	$w = $1;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;


#print '%', ("123456789^" x 10), "\n";
my %ref;
my $hits;
while(<>) {
    #print "$_";
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit);
    # Header lines contain no hits
    unless (%hit) {
	# Print header to maintain a valid SAM output
#	print "$_\n";
	next;
    }
    # Nothing to do if there is no hit for the query sequence, so skip it
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
    &plothit(\%hit, \%ref, $w);
    $hits++;
}
print "No hits to plot in the SAM file\n" unless $hits;

sub plothit{
    # Prints an ASCII plot for the hit to STDOUT
    my ($href, $ref, $width, $flipped) = @_;
    my %hit = %$href;
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});

    # Get R info
    my $start = $hit{'POS'};
    my $aligned = $aln->{'length'} - $aln->{'insertion'};
    my $len = $ref->{ $hit{'RNAME'} };

    # Get Q info
    my $q = $aln->{'length'} - $aln->{'deletion'};
    my $qlen = $q + $aln->{'unmapped'};

    # Calculate postions
    my $a;
    my $b;
    $a = sprintf("%.0f", ($start -1) / $len * $width);
    $b = sprintf("%.0f", $aligned / $len * $width);
    # Add a point even if too small
    $b = 1 unless $b;
    my $c = $width - ($a + $b);
    if ($c < 0) {
	# If 45.5 is gap (a) and 54.5 is match (b) then it will be rounded
	#  to 46 and 55. And that is longer than the width.
	# Adjsut so it still fits and give preference to match over gap.
	$a--;
	$c++;
    }

    # Calculate identity
    my $identity = ($aln->{'length'} - $hit{'NM:i'}) / $aln->{'length'};
    
    # Plot
    my $match = ">";
    $match = "<" if $hit{'FLAG'} & 16;
    my $gap = "_";
    print "|" . ($gap x $a) . ($match x $b) . ($gap x $c) . "|";
    my $qblock = $hit{'QNAME'} . " (". sprintf("%.2f", $q / $qlen * 100) ."\%; l=" . $qlen . ")";
    my $rblock = $hit{'RNAME'} . " (". sprintf("%.2f", $aligned / $len * 100) ."\%; l=" . $len . ")";
    if ($flipped) {
	print " R:$qblock Q:$rblock";
    } else {
	print " R:$rblock Q:$qblock";
    }
    print " I:" . sprintf("%.4f", $identity) ."\n";
}
