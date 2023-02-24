#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointbasics;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool convert to SAM format files to (nucmer) delta format printed on STDOUT.\n" . 
    "\tDisclaimer: the output lacks the header line that specifies the input files\n" .
    "\tfor the seqence comparison (in this case the two files that contain the \n" .
    "\treference sequences and the query sequences respectively).\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] [SAM file]\n";
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

my $ref;
my $query;

my %reflen;

my $r;
my $t;
my @mem;
my $cigar;
my $bitflag;
my $pos;
my $length;
my $covered;
my $tlen;
my $edit;
my $seq = "";
my $qual = "";
my $mapq;
my $x;
my $y;
my @rest;
my $pair = {};
while(<>) {
    if (/^\@/) {
	if (/^\@SQ\s+SN:(\S+)\s+LN:(\d+)/) {
	    $reflen{$1} = $2;
	} elsif (/^\@PG\s+ID:(\S+)\s+PN:(\S+)\s+VN:(\S+)\s+CL:(.+)/) {
	    my ($id, $pn, $vn, $cl) = ($1, $2, $3, $4);
	    # Not sure what to do with these
	}
	next;
    }
    s/\R//g;
    my @tab = split/\t/;
    ($t, $bitflag, $r, $pos, $mapq, $cigar, $x, $y, $tlen, $seq, $qual, @rest) = @tab;
    /NM:i:(\d+)/;
    my $edit = $1;
    
    # lengths are missing
    
    my $mem = 0;
    my $tlen = 0;
    my $delta = "";
    my $rcov = 0;
    my $tcov = 0;
    my $t1 = 1;
    # Remove clipped ends and adjust length and start position
    if ($cigar =~ s/^(\d+)[SH]//) {
	$t1 += $1 unless $bitflag & 16;
	$tlen += $1;
    }
    if ($cigar =~ s/(\d+)[SH]$//) {
	$t1 += $1 if $bitflag & 16;
	$tlen += $1;
    }
    while (length $cigar) {
	$cigar =~ s/^(\d+)//;
	my $length = $1;
	$cigar =~ s/^(\D)//;
	my $type = $1;
	if ($type eq "D" || $type eq "N") {
	    $rcov += $length;
	    $delta .= ($mem + 1) . "\n";
	    if ($length > 1) {
		$delta .= "1\n" x ($length - 1);
	    }
	} elsif ($type eq "I") {
	    $tcov += $length;
	    $tlen += $length;
	    $delta .= "-" . ($mem + 1) . "\n";
	    if ($length > 1) {
		$delta .= "-1\n" x ($length - 1);
	    }
	} elsif ($type eq "M" || $type eq "=" || $type eq "X") {
	    $rcov += $length;
	    $tcov += $length;
	    $tlen += $length;
	    $mem = $length;
	} else {
	    print STDERR "Unexpected cigar element found ('$length$type') that will be skipped\n";
	}
    }
    $delta .= "0\n";
    # matched reagions needed
    $pair->{$r}->{$t}->{'head'} = ">$r $t $reflen{$r} $tlen\n";
    # print ">$r $t $reflen{$r} $tlen\n";
    my $t2 = $t1 + $tcov - 1 ;
    if ($bitflag & 16) {
	($t2, $t1) = ($t1, $t2);
    }
    $pair->{$r}->{$t}->{'body'} .= join(" ", $pos, $pos + $rcov - 1, $t1, $t2, $edit, $edit, "0\n") . $delta;
    # print join(" ", $pos, $pos + $rcov - 1, $t1, $t2, $edit, $edit, "0\n");
    # print $delta;
}
for my $r (sort keys %$pair) {
    for my $t (sort keys %{ $pair->{$r} }) {
	print $pair->{$r}->{$t}->{'head'};
	print $pair->{$r}->{$t}->{'body'};
    }
}
