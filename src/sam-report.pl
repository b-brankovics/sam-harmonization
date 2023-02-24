#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use Data::Dumper;
use biointbasics;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n" .
    "\tA tool to summarize all the hits from a SAM file. The tool either opens the file\n" .
    "\tspecified as input or reads from STDIN when no file is given. The program prints a\n" .
    "\tpairwise summary for refernce-query pairs to STDOUT as tab-separated data.\n";
my $usage = 
    "Usage:\n\t$0 [-h | --help] [SAM file]";
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


my $filter;
my @keep;
for (@ARGV) {
    if (/^--?b(est)?$/) {
	$filter++;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;


my %ref;
my %query;
my %lengths;

# Keep all the individual identity scores
my @scores;

# Store pairwise information by collecting all hits between two IDs
my $r2q = {};
my $q2r = {};

while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    next unless %hit;
    if ($hit{'RNAME'} eq "*") {
	$query{ $hit{'QNAME'} } = length($hit{'SEQ'});
	next;
    }
    # Get alginment data
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    # Store the length of queey seq
    $query{ $hit{'QNAME'} } = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};

    my $i = $hit{'RNAME'};
    my $j = $hit{'QNAME'};
    # Store pair information
    #  Link if already exists
    #  Otherwise create it
    my $pair;
    if ($r2q->{$i}->{$j}) {
	$pair = $r2q->{$i}->{$j};
    } else {
	$pair = {};
	$r2q->{$i}->{$j} = $pair;
	$q2r->{$j}->{$i} = $pair;
	$pair->{"A"}->{"hits"} = {};
	$pair->{"B"}->{"hits"} = {};
	
    }

    my $score = ($aln->{'length'} - $hit{'NM:i'}) / $aln->{'length'};
    $pair->{"identical"} += $aln->{'length'} - $hit{'NM:i'};
    $pair->{"diff"} += $hit{'NM:i'};
    $pair->{"aligned"} += $aln->{'length'};
    push @{ $pair->{"hits"} }, \%hit;
    push @{ $pair->{"scores"} }, $score;

    # store info for coverage
    # ref
    my $s = $hit{'POS'};
    my $e = $hit{'POS'} + $aln->{'length'} - $aln->{'insertion'} - 1;
    my $r = $pair->{"A"}->{"hits"};
    if ($r->{$i} && $r->{$i}->{$s}) {
	if ($r->{$i}->{$s} < $e) {
	    $r->{$i}->{$s} = $e;
	}
    } else {
	$r->{$i}->{$s} = $e;
    }
    #  for query
    my $s2 = $aln->{'start'};
    my $e2 = $aln->{'start'} - 1 + $aln->{'length'} - $aln->{'deletion'};
    my $q = $pair->{"B"}->{"hits"};
    if ($q->{$j} && $q->{$j}->{$s2}) {
	if ($q->{$j}->{$s2} < $e2) {
	    $q->{$j}->{$s2} = $e2;
	}
    } else {
	$q->{$j}->{$s2} = $e2;
    }
}

#print "# ", join("\t", "A-ID", "B-ID", "Similarity", "Identical", "Aligned", "A-coverage", "B-coverage", "A-covered", "B-covered", "A-length", "B-length", "orientation", "A-region", "B-region"), "\n";
print "# ", join("\t", "A-ID", "B-ID", "Similarity", "Identical", "Aligned", "A-coverage", "B-coverage", "A-covered", "B-covered", "A-length", "B-length", "A-region", "B-region"), "\n";
# Preprocess pairs
for my $i (sort keys %ref) {
    next unless $r2q->{$i};
    for my $j (sort keys %query) {
	next unless $r2q->{$i}->{$j};
	my $pair = $r2q->{$i}->{$j};
	$pair->{"A"}->{"region"} = [];
	$pair->{"B"}->{"region"} = [];
	$pair->{"A"}->{"covered"} = &coverage($pair->{"A"}->{"hits"}, $pair->{"A"}->{"region"});
	$pair->{"B"}->{"covered"} = &coverage($pair->{"B"}->{"hits"}, $pair->{"B"}->{"region"});
	$pair->{"similarity"} = $pair->{"identical"} / $pair->{"aligned"};
	$pair->{"max"} = $pair->{"identical"} / (1 + $pair->{"diff"});
	# my $cov_a = $cov1 / $ref{$i};
	# my $cov_b = $cov2 / $query{$j};

	# print Dumper($pair->{"B"}->{"hits"}), "\n";
	# print Dumper($pair->{"B"}->{"region"}), "\n";
	# print Dumper($pair->{"A"}->{"covered"}), "\n";
    }
}	

my %found;

# Best for each query
my @best;
for my $j (sort keys %query) {
    next unless $q2r->{$j};
#    my $term = "similarity";
    my $term = "max";
    my ($i) = sort{$q2r->{$j}->{$b}->{$term} <=> $q2r->{$j}->{$a}->{$term}} keys %{ $q2r->{$j} };
    push @best, $q2r->{$j}->{$i};
    $q2r->{$j}->{$i}->{"A"}->{"ID"} = $i;
    $q2r->{$j}->{$i}->{"B"}->{"ID"} = $j;
}
if ($filter) {
    for my $pair (@best) {
	my $i = $pair->{"A"}->{"ID"};
	my $j = $pair->{"B"}->{"ID"};
	$found{$i}++;
	$found{$j}++;
	    print join("\t",
		       $i,
		       $j,
		       $pair->{"identical"} / $pair->{"aligned"},
		       $pair->{"identical"},
		       $pair->{"aligned"},
		       $pair->{"A"}->{"covered"} / $ref{$i},
		       $pair->{"B"}->{"covered"} / $query{$j},
		       $pair->{"A"}->{"covered"},
		       $pair->{"B"}->{"covered"},
		       $ref{$i},
		       $query{$j},
		       # $ori,
		       join(";", @{ $pair->{"A"}->{"region"} }),
		       join(";", @{ $pair->{"B"}->{"region"} })
		), "\n";
    }
} else {
    for my $i (sort keys %ref) {
	next unless $r2q->{$i};
	for my $j (sort keys %query) {
	    next unless $r2q->{$i}->{$j};
	    my $pair = $r2q->{$i}->{$j};
	    $found{$i}++;
	    $found{$j}++;
	    print join("\t",
		       $i,
		       $j,
		       $pair->{"identical"} / $pair->{"aligned"},
		       $pair->{"identical"},
		       $pair->{"aligned"},
		       $pair->{"A"}->{"covered"} / $ref{$i},
		       $pair->{"B"}->{"covered"} / $query{$j},
		       $pair->{"A"}->{"covered"},
		       $pair->{"B"}->{"covered"},
		       $ref{$i},
		       $query{$j},
		       # $ori,
		       join(";", @{ $pair->{"A"}->{"region"} }),
		       join(";", @{ $pair->{"B"}->{"region"} }),
		       $pair->{"max"}
		), "\n";
	    
	}
    }
}
# Report the ones that had no hits
for my $i (sort keys %ref) {
    next if $found{$i};
    print join("\t",
	       $i,
	       "NA",
	       "NA", #$pair->{"identical"} / $pair->{"aligned"},
	       "NA", #$pair->{"identical"},
	       "NA", #$pair->{"aligned"},
	       0, #$pair->{"A"}->{"covered"} / $pair->{"A"}->{"length"},
	       "NA", #$pair->{"B"}->{"covered"} / $pair->{"B"}->{"length"},
	       0, #$pair->{"A"}->{"covered"},
	       "NA", #$pair->{"B"}->{"covered"},
	       $ref{$i},
	       "NA",
	       # "NA",
	       "NA",
	       "NA"
	), "\n";
}    
for my $j (sort keys %query) {
    next if $found{$j};
    print join("\t",
	       "NA",
	       $j,
	       "NA", #$pair->{"identical"} / $pair->{"aligned"},
	       "NA", #$pair->{"identical"},
	       "NA", #$pair->{"aligned"},
	       "NA", #$pair->{"A"}->{"covered"} / $pair->{"A"}->{"length"},
	       0, #$pair->{"B"}->{"covered"} / $pair->{"B"}->{"length"},
	       "NA", #$pair->{"A"}->{"covered"},
	       0, #$pair->{"B"}->{"covered"},
	       "NA",
	       $query{$j},
	       # "NA",
	       "NA",
	       "NA"
	), "\n";
}    
    



# subroutine
sub coverage {
    my ($r, $region) = @_;
    # Merge overlaps or continuous stretches for reference
    for my $id (sort keys %$r) {
	my $from;
	my $to;
	for my $s (sort{$a <=> $b} keys %{ $r->{$id} }) {
	    if ($to && $to + 1 >= $s) {
		if ($r->{$id}->{$s} > $to) {
		    # extend the old one
		    $r->{$id}->{$from} = $r->{$id}->{$s};
		}
		# remove the current one
		delete $r->{$id}->{$s};
		$s = $from;
	    }
	    $from = $s;
	    $to = $r->{$id}->{$s};
	}
    }

    # Sum total coverage of reference
    my $covr;
    for my $id (sort keys %$r) {
	my @range;
	for my $s (sort{$a <=> $b} keys %{ $r->{$id} }) {
	    $covr += $r->{$id}->{$s} + 1 - $s;
	    push @range, $s . ".." . $r->{$id}->{$s};
	}
	if (@range) {
	    my $location = $id . ":";
	    $location .= "join(" if 1 < scalar @range;
	    $location .= join(",", @range);
	    $location .= ")" if 1 < scalar @range;
	    push @$region, $location;
	}
    }
    return $covr;
}


