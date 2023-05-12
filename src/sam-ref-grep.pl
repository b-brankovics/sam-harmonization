#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

my $programname = "sam-ref-grep.pl";
my $version = "1.0";
my $cmd = join(" ", $programname, @ARGV);

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to filter SAM files using grep like filter based on reference names.\n" .
    "\tThe tool either opens the file specified as input or reads from STDIN when no file is given.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] [SAM file] ['<regex>']\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-v | --invert-match\n\t\tInvert the sense of matching, to select non-matching entries.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

my %ref;
my %query;

my $invert;
my @requests;
my @keep;
my $input;

# Print help if needed
biointbasics::print_help(\@ARGV, $info);


for (@ARGV) {
    if (/^-$/ || /\.sam$/) {
	push @keep, $_;
    } elsif (/^-v$/) {
	$invert = "yes";
#	print STDERR "Using inverted mode: only those lines are reported that do not pass\n";
    } else {
	push @requests, $_;
    }
}
@ARGV = @keep;

biointbasics::print_help(\@ARGV, $info, "ERROR: No regular expressions specified for filtering\n\n") unless @requests;

my $header = "true";
my $previousprogram = "";
while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    unless (%hit) {
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

    my $positive;
#    $positive++ if $filter{ $hit{'RNAME'} };
    for my $regex (@requests) {
	$positive++ if $hit{'RNAME'} =~ /$regex/;
    }
        
    if (($positive && ! $invert) || (! $positive && $invert)) {
	# Print if it matched and not inverted mode
	# Print if it did not match but invert is selected
	print biointsam::sam_string(\%hit), "\n";
    }
}
