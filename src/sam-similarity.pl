#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;


#===DESCRIPTION=================================================================
# A tool to remove duplicate sequence from FASTA format files and
#  print the groups to STDERR

my $description = 
    "Description:\n\t" .
    "A tool to calculate ANI and coverage statistics for SAM files.\n" .
    "\tThe tool either opens the file specified as input or reads from STDIN when no file is given.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS]  [SAM file]\n";
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

biointbasics::print_help(\@ARGV, $info);

my %ref;
my %query;

my $files = 1;

# ANI calculation
my $numer; # nummerator
my $denom; # denominator
my $length;
my $mis;
# Store homology stretches per input
my %r;
my $i;
my %q;
my $j;
# Keep all the individual identity scores
my @scores;


while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    next unless %hit;
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
    # Get alginment data
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    # Store the length of queey seq
    $query{ $hit{'QNAME'} } = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};

    $i = $hit{'RNAME'};
    $j = $hit{'QNAME'};
    my $length = $aln->{'length'};
    my $mis = $hit{'NM:i'};
    $numer += $length - $mis;
    $denom += $length;
    push @scores, ($length - $mis) / $length;
    # store info for coverage
    # ref
    my $s = $hit{'POS'};
    my $e = $hit{'POS'} + $length - $aln->{'insertion'} - 1;
    if ($r{$i} && $r{$i}->{$s}) {
	if ($r{$i}->{$s} < $e) {
	    $r{$i}->{$s} = $e;
	}
    } else {
	$r{$i}->{$s} = $e;
    }
    # q
    my $s2 = $aln->{'start'};
    my $e2 = $aln->{'start'} - 1 + $aln->{'length'} - $aln->{'deletion'};
    if ($q{$j} && $q{$j}->{$s2}) {
	if ($q{$j}->{$s2} < $e2) {
	    $q{$j}->{$s2} = $e2;
	}
    } else {
	$q{$j}->{$s2} = $e2;
    }
}



# Merge overlaps or continuous stretches for reference
for my $id (sort keys %r) {
    my $from;
    my $to;
    for my $s (sort{$a <=> $b} keys %{ $r{$id} }) {
	if ($to && $to + 1 >= $s) {
	    if ($r{$id}->{$s} > $to) {
		# extend the old one
		$r{$id}->{$from} = $r{$id}->{$s};
	    }
	    # remove the current one
	    delete $r{$id}->{$s};
	    $s = $from;
	}
	$from = $s;
	$to = $r{$id}->{$s};
    }
}

# Sum total coverage of reference
my $covr;
for my $id (sort keys %r) {
    for my $s (sort{$a <=> $b} keys %{ $r{$id} }) {
	$covr += $r{$id}->{$s} + 1 - $s;
    }
}

# Merge overlaps or continuous stretches for query
for my $id (sort keys %q) {
    my $from;
    my $to;
    for my $s (sort{$a <=> $b} keys %{ $q{$id} }) {
	if ($to && $to + 1 >= $s) {
	    if ($q{$id}->{$s} > $to) {
		# extend the old one
		$q{$id}->{$from} = $q{$id}->{$s};
	    }
	    delete $q{$id}->{$s};
	    $s = $from;
	}
	$from = $s;
	$to = $q{$id}->{$s};
    }
}

# Sum total coverage of query
my $covq;
for my $id (sort keys %q) {
    for my $s (sort{$a <=> $b} keys %{ $q{$id} }) {
	$covq += $q{$id}->{$s} + 1 - $s;
    }
}

# Report ANI and similarity
unless ($denom) {
    print "No alignment in the SAM file\n";
    exit;
}

print "ANI: " . (sprintf '%.2f', 100 * $numer / $denom) . "% ($numer/$denom)\n";
 


my @sort = sort {$a<=>$b} @scores;
my $sum;
for (@sort) {
    $sum += $_;
}
print "(Range of similarity scores [$sort[0], $sort[-1]]; arithmetic mean " . ($sum / scalar @sort) . ")\n";
# print "alingment length: $denom\n";
# print "identical positions: $numer\n";

# Get full lengths and report coverage info
my ($len_r, $len_q);
for (values %ref) {
    $len_r += $_;
}
for (values %query) {
    $len_q += $_;
}

if ($files) {
    print "R covered: $covr";
    print " (", (sprintf '%.2f', (100 * $covr / $len_r)), "%)";
    print "\n";
    print "Q covered: $covq";
    print " (", (sprintf '%.2f', (100 * $covq / $len_q)), "%)";
    print "\n";
} else {
    $len_r = ">=$len_r";
    $len_q = ">=$len_q";
}


print "# ", join("\t", "Similarity", "Identical", "Aligned", "A-covered", "B-covered", "A-length", "B-length"), "\n";
print join("\t", $numer / $denom, $numer, $denom, $covr, $covq, $len_r, $len_q), "\n";


