#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use Data::Dumper;

#===DESCRIPTION=================================================================
# A tool to remove duplicate sequence from FASTA format files and
#  print the groups to STDERR

my $description = 
    "Description:\n" .
    "\tA tool to summarize all the hits from a SAM file. The tool either opens the file\n" .
    "\tspecified as input or reads from STDIN when no file is given. The program produces\n" .
    "\ttwo different outputs: it prints a pairwise summary for refernce-query pairs to STDOUT\n" .
    "\tas tab-separated data, and it also prints summary per reference entry to STDERR, also as TSV.\n";
my $usage = 
    "Usage:\n\t$0 [-h | --help] [-b | --best]  [SAM file] >pairwise-summary.tsv 2>per-reference-summary.tsv\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-b | --best\n\t\tKeep only one reference for each query sequence.\n" .
    "\n";

#===MAIN========================================================================

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

# Print help if needed
&print_help(\@ARGV);

my %ref;
my %query;
my %lengths;

# Keep all the individual identity scores
my @scores;

# Store pairwise information by collecting all hits between two IDs
my $r2q = {};
my $q2r = {};

# Coverage depth
# depth -> RefID -> start -> end -> coverage depth
my $depth = {};

while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    next unless %hit;
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
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
    $pair->{"rawscore"} += $hit{'AS:i'} if $hit{'AS:i'};
    push @{ $pair->{"hits"} }, \%hit;
    push @{ $pair->{"scores"} }, $score;

    # store info for coverage (breadth)
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
    # For ref coverage depth
    &cov_depth($i, $s, $e, $depth);
    
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


my @list = sort keys %ref;
my $ref_r2q = $r2q;

if ($filter) {
    # Best for each query
    my @best;
    my $keep = {};
    for my $j (sort keys %query) {
	next unless $q2r->{$j};
	# my $term = "similarity";
	# my $term = "max";
	my $term = "rawscore";
	# Select the best
	my ($i) = sort{$q2r->{$j}->{$b}->{$term} <=> $q2r->{$j}->{$a}->{$term}} keys %{ $q2r->{$j} };
	push @best, $q2r->{$j}->{$i};
	$keep->{$i}->{$j} = $q2r->{$j}->{$i};
	$q2r->{$j}->{$i}->{"A"}->{"ID"} = $i;
	$q2r->{$j}->{$i}->{"B"}->{"ID"} = $j;
    }
    @list = @best;
    $ref_r2q = $keep;
}

# Summary for each ref
#  If best match filtering was used then hash reference is already corrected for it
print STDERR join("\t", "# RefID", "ANI", "coverage", "covered", "length", "identical", "aligned", "average depth", "multi covered", "multi coverage", "hits", "total query length", "total query mapped"), "\n";
for my $i (sort keys %ref) {
    # Sum up query info: count, total length, total mapped regions length
    my $count = 0;
    my $total = 0;
    my $mapped = 0;
    my $sum_covered = 0;
    if ($ref_r2q->{$i}) {
	# collect all regions
	# collect identical
	# collect aligned lengths
	my $r = {};
	my $identical = 0;
	my $length = 0;
	for my $j (sort keys %{ $ref_r2q->{$i} }) {
	    $count++;
	    my $pair = $ref_r2q->{$i}->{$j};
	    $identical += $pair->{"identical"};
	    $length += $pair->{"aligned"};
	    $total += $query{$j};
	    $mapped += &coverage($pair->{"B"}->{"hits"}, []);
	    # store info for coverage
	    for my $s (keys %{ $pair->{"A"}->{"hits"}->{$i} }) {
		my $e = $pair->{"A"}->{"hits"}->{$i}->{$s};
		# Add length of ref covered by alignment
		$sum_covered += $e - $s + 1;
		if ($r->{$i} && $r->{$i}->{$s}) {
		    if ($r->{$i}->{$s} < $e) {
			$r->{$i}->{$s} = $e;
		    }
		} else {
		    $r->{$i}->{$s} = $e;
		}
	    }
	}
	# Calculate values
	my $regions = [];
	my $covered = &coverage($r, $regions);
	my $ani = $identical / $length;
	my $avg_depth = $sum_covered / $covered;
	my $multicov = 0;
	for my $s (keys %{ $depth->{$i} }) {
	    my ($e) = keys %{ $depth->{$i}->{$s} };
	    $multicov += $e - $s + 1 if $depth->{$i}->{$s}->{$e} > 1;
	}
	print STDERR join("\t", $i, $ani, sprintf("%.4f", ( $covered / $ref{$i} )), $covered, $ref{$i}, $identical, $length, $avg_depth, $multicov, $multicov / $covered, $count, $total, $mapped), "\n";
	#print STDERR join("\t", $i, $ani, sprintf("%.4f", ( $covered / $ref{$i} )), $covered, $ref{$i}, $identical, $length, $avg_depth, $count, $total, $mapped), "\n";
    } else {
	print STDERR join("\t", $i, "NA", "NA", "NA", $ref{$i}, "NA", "NA", "NA", "NA", "NA", $count, "NA", "NA"), "\n";
    }
}

# Print info for reference and query pairs
#  If best match filtering was used then hash reference is already corrected for it
# Print header
print "# ", join("\t", "A-ID", "B-ID", "Similarity", "Identical", "Aligned", "A-coverage", "B-coverage", "A-covered", "B-covered", "A-length", "B-length", "A-region", "B-region"), "\n";
# Mark sequences that have hits, so later we can check what needs to be printed as NA
my $found = {"r" => {}, "q" => {}};

for my $i (sort keys %$ref_r2q) {
    next unless $ref_r2q->{$i};
    for my $j (sort keys %query) {
	next unless $ref_r2q->{$i}->{$j};
	my $pair = $ref_r2q->{$i}->{$j};
	#    my $i = $pair->{"A"}->{"ID"};
	#    my $j = $pair->{"B"}->{"ID"};
	$found->{"r"}->{$i}++;
	$found->{"q"}->{$j}++;
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
}
# Report the ones that had no hits
for my $i (sort keys %ref) {
    next if $found->{"r"}->{$i};
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
    next if $found->{"q"}->{$j};
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


sub cov_depth {
    my ($i, $s, $e, $depth) = @_;
    # For ref coverage depth
    # check overlap
    # clip overlap and keep rest if any
    for my $old_start (sort{$a<=>$b} keys %{ $depth->{$i} }) {
	# Stop of all used up
	last unless ($s <= $e);
	# There can only be one end for each region
	my ($old_end) = keys %{ $depth->{$i}->{$old_start} };
	my $old_cov = $depth->{$i}->{$old_start}->{$old_end};
	if ($old_start < $s) {
	    if ($s <= $old_end) {
		if ($e < $old_end) { # nested, old splits into 3
		    # Clip old: create adjusted range and delete old one
		    $depth->{$i}->{$old_start}->{$s - 1} = $old_cov;
		    $depth->{$i}->{$e + 1}->{$old_end} = $old_cov;
		    delete $depth->{$i}->{$old_start}->{$old_end};
		    # Create new for overlap region
		    $depth->{$i}->{$s}->{$e} = $old_cov + 1;
		    # current used up
		    $s = $e + 1;
		} elsif ($e == $old_end) { # new covered, old splits in 2
		    # Clip old: create adjusted range and delete old one
		    $depth->{$i}->{$old_start}->{$s - 1} = $old_cov;
		    delete $depth->{$i}->{$old_start}->{$old_end};
		    # Create new for overlap region
		    $depth->{$i}->{$s}->{$e} = $old_cov + 1;
		    # current used up
		    $s = $e + 1;
		} else { # old, overlap, new
		    # examine and add and clip
		    my $overlap = $old_end - $s + 1;
		    # Clip old: create adjusted range and delete old one
		    $depth->{$i}->{$old_start}->{$old_end - $overlap} = $old_cov;
		    delete $depth->{$i}->{$old_start}->{$old_end};
		    # Create new for overlap region
		    $depth->{$i}->{$old_end - $overlap + 1}->{$old_end} = $old_cov + 1;
		    # Modify current
		    $s += $overlap;
		}
	    } 
	} elsif ($old_start == $s) {
	    if ($e < $old_end) {
		my $overlap = $e - $s + 1;
		# Clip old: create adjusted range and delete old one
		$depth->{$i}->{$old_start}->{$old_start + $overlap - 1} = $old_cov + 1;
		delete $depth->{$i}->{$old_start}->{$old_end};
		# Create new for overlap region
		$depth->{$i}->{$old_start + $overlap}->{$old_end} = $old_cov;
		# current is used up
		$s = $e + 1;
	    } elsif ($e == $old_end) {
		$depth->{$i}->{$old_start}->{$old_end} += 1;
		# current is used up
		$s = $e + 1;
	    } else {
		# $old_end < $e
                # Old one is completely covered
                $depth->{$i}->{$old_start}->{$old_end} = $old_cov + 1;
                # Modify current
		$s = $old_end + 1;
	    }
	} elsif ($old_start <= $e) {
	    # situation $s < $old_start <= $e
	    if ($old_end < $e) {
		# old is completely covered
		$depth->{$i}->{$old_start}->{$old_end} = $old_cov + 1;
		# add a new before the old
		$depth->{$i}->{$s}->{$old_start - 1} = 1;
		# adjust current
		$s = $old_end + 1;
	    } elsif ($old_end == $e) {
		# old is completely covered
		$depth->{$i}->{$old_start}->{$old_end} = $old_cov + 1;
		# adjust current
		$e = $old_start - 1;
	    } else {
		# situation $s < $old_start <= $e < $old_end
		my $overlap = $e - $old_start + 1;
		# Clip old: create adjusted range and delete old one
		$depth->{$i}->{$old_start}->{$e} = $old_cov + 1;
		delete $depth->{$i}->{$old_start}->{$old_end};
		# Create new region
		$depth->{$i}->{$e + 1}->{$old_end} = $old_cov;
		# adjust current
		$e = $old_start - 1;
	    }
	}
    }
    # add the rest
    if ($s <= $e) {
	# Add as region with cov=1
	$depth->{$i}->{$s}->{$e} = 1;
    }
}

