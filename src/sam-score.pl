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
    "A tool to recalculate alignment score for SAM files.\n" .
    "\tThis is useful if hits generated using different tools are being compared.\n" .
    "\tAdditional statistics can also be calculated for hits when needed.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS]  [SAM file]\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-pid | --percent-identity\n\t\tCalculate percent identity (100 * identical positions/alignemnt length).\n" .
    "\t-rid | --reference-identity\n\t\tCalculate reference identity (identical positions/reference length).\n" .
    "\t-qid | --query-identity\n\t\tCalculate query identity (identical positions/query length).\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================


my @keep;
my %add;
for (@ARGV) {
    if (/^--?rid$/ || /^--?reference-identity$/) {
	$add{'RI:f'}++;
    } elsif (/^--?qid$/ || /^--?query-identity$/) {
	$add{'QI:f'}++;
    } elsif (/^--?pid$/ || /^--?percent-identity$/) {
	$add{'PI:f'}++;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;

# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my %ref;

while(<>) {
    # Print header as it is
    print $_ if /^\@/;
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    next unless %hit;
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
    biointsam::score(\%hit);
    if ($add{'RI:f'} || $add{'QI:f'} || $add{'PI:f'}) {
	# Get alginment data
	my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
	# Get identical positions
	my $identical = $aln->{'length'} - $hit{'NM:i'};
	if ($add{'PI:f'}) {
	    $hit{'PI:f'} = sprintf ("%.2f", $identical / $aln->{'length'} * 100);
	}
	if ($add{'RI:f'}) {
	    $hit{'RI:f'} = sprintf ("%.4f", $identical / $ref{ $hit{'RNAME'} });
	}
	if ($add{'QI:f'}) {
	    my $query_len = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};
	    $hit{'QI:f'} = sprintf ("%.4f", $identical / $query_len);
	}
    }
    print biointsam::sam_string(\%hit), "\n";

}
