#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the parent directory
use biointsam;

use File::Basename;
use Data::Dumper;

my $programname = "sam-keep-best.pl";
my $version = "1.2";
my $cmd = join(" ", $programname, @ARGV);

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to filter SAM files and keep the best match/mapping for each query entry and\n" .
    "\tprints to STDOUT. The best match is identified (by default) by summing all the\n" .
    "\talignment scores for all matches belonging to the same reference-query pair.\n" .
    "\tThe tool either opens the file specified as input or reads from STDIN when no file is given.\n";
my $usage = 
    "Usage:\n\t$0 [Options]  [SAM file]\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-g | --global\n\t\tThe best match is identified by summing all the alignment scores for all matches belonging to the same reference-query pair (This is the default mode).\n" .
    "\t-l | --local\n\t\tKeep the best local hits for each region of query sequences.\n\t\tEliminates hits that overlap with another that has a higher score.\n" .
    "\t-mo=<int> | --max-overlap=<int>\n\t\tAllow <int> overlap between hits for local mode. (Default is 0 bp. Incompatible with '-po=<float>'.)\n" .
    "\t-po=<float> | --max-percent-overlap=<float>\n\t\tAllow <float> (10% should be set as '-po=0.1') overlap between hits for local mode. (Incompatible with '-mo=<int>'. Not set by default.)\n" . 
                                                   "\t\tThe overlap has to be shorter than the selected percentage of the shorter hit (understood as matched region on query).\n" .
    "\t-bs | --bitscore\n\t\tInstead of using the alignment score, use the bitscore.\n\t\t(Only possible if the SAM was generated based on a BLAST run)\n" .
    "\t-pi | --percent-identity\n\t\tInstead of using the alignment score, use the percent identity (identical bases/alignment length).\n" .
    "\t-rid | --reference-identity\n\t\tInstead of using the alignment score, use reference identity (identical positions/reference length).\n\t\t(Only implemented for local mode.)\n" .
    "\t-qid | --query-identity\n\t\tInstead of using the alignment score, use query identity (identical positions/query length).\n\t\t(Only implemented for local mode.)\n" .
    "\t-u | --unique\n\t\tKeep only one hit per reference-query pair (or the equally best ones).\n\t\tThis happens prior the other ranking based filtering.\n" .
    "\t-top=<int> | --top=<int>\n\t\tKeep the hits with the <int> highest score. Local mode can only be used with top=1 (the default setting).\n" .
    "\t-n=<int> | --n=<int>\n\t\tPrint only the first <int> hits with the selected high score range (defined by top).\n" .
                             "\t\tBy default n is not set, and all hits with the highest score are printed.\n" .
    "\n";

#===MAIN========================================================================

# Number of top scores to keep. Default is 1.
#   Cannot be higher than 1 for local mode
my $top = 1;
# How many to keep in total?
#   - default will keep all that fit to the top selection
#   - if set then only report the "first" n hits in the ranking
my $n;
# Modes:
#   - summative (global)
#   - local
my $mode = "summative";
# Ranking criterion:
#   - AS:i (alignment score)
#   - BS:i (bit score)
#   - PI:f (percent identity/similarity: identical/alignment length)
#   - RI:f (similarity to the reference sequences)
#   - QI:f (similarity to the query sequences)
my $unique;
my $term = 'AS:i';
# Maximal overlap (for local mode):
#   default is 0 bp
my $maxoverlap = 0;
my $percentoverlap = 0;

my @keep;
for (@ARGV) {
    if (/^--?bs$/ || /^--?bitscore$/) {
	$term = 'BS:f';
	print STDERR "Switched to using bitscore ('BS:f') instead of alignment score ('AS:i').\n";
    } elsif (/^--?g$/ || /^--?global$/) {
	$mode = "summative";
    } elsif (/^--?l$/ || /^--?local$/) {
	$mode = "local";
    } elsif (/^--?u$/ || /^--?unique$/) {
	$unique++;
    } elsif (/^--?pi$/ || /^--?(percent-)?identity$/) {
	$term = 'PI:f';
	print STDERR "Switched to using percent identity ('PI:f') instead of alignment score ('AS:i').\n";
    } elsif (/^--?rid$/ || /^--?reference-identity$/) {
	$term = 'RI:f';
	$mode = "local";
	print STDERR "Switched to using reference identity ('RI:f') instead of alignment score ('AS:i').\n";
	print STDERR "Switched to using local mode instead of summative.\n";
    } elsif (/^--?qid$/ || /^--?query-identity$/) {
	$term = 'QI:f';
	$mode = "local";
	print STDERR "Switched to using query identity ('QI:f') instead of alignment score ('AS:i').\n";
	print STDERR "Switched to using local mode instead of summative.\n";
    } elsif (/^--?top=(\d+)$/) {
	$top = $1;
	print STDERR "Keeping the top $top hits instead of only the ones with the highest score.\n" if $top > 1;
    } elsif (/^--?n=(\d+)$/) {
	$n = $1;
	print STDERR "Printing the first $n results for the top $top hits instead of all of those. (Returns 0-$n hits per case.)\n";
    } elsif (/^--?mo=(\d+)$/ || /^--?max-overlap=(\d+)$/) {
	$maxoverlap = $1;
	print STDERR "Overlap tolerance is increased to $maxoverlap bp from 0.\n";
    } elsif (/^--?po=(\d\.\d+)$/ || /^--?max-percent-overlap=(\d+)$/) {
	$percentoverlap = $1;
	print STDERR "Overlap tolerance is increased to " . sprintf("%0.2f", (100 * $percentoverlap)) ." % percent (of the shorter overlaping segment) from 0 bp.\n";
    } elsif ($_ eq '-') {
	push @keep, $_;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;

# Print help if needed
&print_help(\@ARGV);

if ($top > 1 && $mode eq "local") {
    &print_help(\@ARGV, "ERROR: Local mode can only be run when 'top' is set to '1' (default setting for 'top').\n");
}
if ($maxoverlap && $mode ne "local") {
    &print_help(\@ARGV, "Max (percent) overlap is only compatible with local mode.\n");
}
if ($percentoverlap && $mode ne "local") {
    &print_help(\@ARGV, "Max (percent) overlap is only compatible with local mode.\n");
}
if ($percentoverlap && $maxoverlap) {
    &print_help(\@ARGV, "Max overlap and max percent overlap cannot be used at the same time.\n");
}


# print STDERR "@ARGV\n";

my %ref;
my %query;
my %lengths;

# Keep all the individual identity scores
my @scores;

# Store pairwise information by collecting all hits between two IDs
my $r2q = {};
my $q2r = {};

# Store information for local best hits
my $location = {};

my $header = "true";
my $previousprogram = "";
while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    unless (%hit) {
	if (/^\@PG\tID:(\S+)/) {
	    $previousprogram = $1;
	}
	print "$_\n";
	next;
    }

    if ($header) { # First line after the header section
	my $text = "\@PG\tID:$programname\tPN:$programname";
	$text .= "\tPP:$previousprogram" if $previousprogram;
	print $text . "\tVN:$version\tCL:$cmd\n";
	$header = undef;
    }

    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";

    # Add read pair flags ('/1' and '/2' for forward and reverse reads of a pair)
    #  otherwise false conclusions will be used
    #  these flags will be removed before printing to OUTPUT
    if ($hit{"FLAG"} & 64) {
	$hit{'QNAME'} .= "/1";
    } elsif ($hit{"FLAG"} & 128) {
	$hit{'QNAME'} .= "/2";
    }

    
    # Get alginment data
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    # Store the length of queey seq
    $query{ $hit{'QNAME'} } = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};

    my $i = $hit{'RNAME'};
    my $j = $hit{'QNAME'};

    # Check (and calculate if needed) identity scores
    if ($term =~ /[PRQ]I:f/ && ! $hit{$term}) {
	my $identical = $aln->{'length'} - $hit{'NM:i'};
	if ($term eq 'QI:f') {
	    $hit{$term} = sprintf ("%.4f", $identical / $query{ $hit{'QNAME'} });
	} elsif ($term eq 'RI:f') {
	    $hit{$term} = sprintf ("%.4f", $identical / $ref{ $hit{'RNAME'} });
	} elsif ($term eq 'PI:f') {
	    $hit{$term} = sprintf ("%.2f", $identical / $aln->{'length'} * 100);
	}
    }
    
    if ($mode eq "local") {
	# Need (IDs, score,) begin and end
	my $from = $aln->{'start'};
	my $to = $aln->{'start'} - 1 + $aln->{'length'} - $aln->{'deletion'};
	my $score = $hit{$term};
	push @{ $location->{$j} }, {
	    from => $from,
	    to => $to,
	    ref => $i,
	    score => $score,
	    hit => \%hit,
	};
    } else {
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
		
	&update_pair($pair, \%hit, $aln, $i, $j, $term);
    }
}


# Filter best unique combos

if ($unique) {
    if ($mode eq "local") {
	for my $i (sort keys %$location) {
	    my %matches;
	    for ( @{ $location->{$i} } ) {
		push @{ $matches{ $_->{'ref'} } }, $_;
	    }
	    my @keep;
	    for my $j (keys %matches) {
		my $max = 0;
		my @best;
		for (@{ $matches{$j} }) {
		    my $score = $_->{'score'};
		    if ($score > $max) {
			@best = ($_);
			$max = $score;
		    } elsif ($score == $max) {
			push @best, $_;
		    }
		}
	 	push @keep, @best;
	    }
	    $location->{$i} = \@keep;
	}
    } elsif ($mode eq "summative") {
	
	#print "# ", join("\t", "A-ID", "B-ID", "Similarity", "Identical", "Aligned", "A-coverage", "B-coverage", "A-covered", "B-covered", "A-length", "B-length", "orientation", "A-region", "B-region"), "\n";
	# Preprocess pairs
	for my $i (sort keys %ref) {
	    next unless $r2q->{$i};
	    for my $j (sort keys %query) {
		next unless $r2q->{$i}->{$j};
		my $pair = $r2q->{$i}->{$j};
		# Select best
		my $max = 0;
		my @best;
		for (@{ $pair->{"hits"} } ) {
		    my $score = $_->{$term};
		    if ($score > $max) {
			@best = ($_);
			$max = $score;
		    } elsif ($score == $max) {
			push @best, $_;
		    }
		}
		#$pair->{"hits"} = \@best;
		# Need to reset some values
		$pair = {};
		$r2q->{$i}->{$j} = $pair;
		$q2r->{$j}->{$i} = $pair;
		$pair->{"A"}->{"hits"} = {};
		$pair->{"B"}->{"hits"} = {};
		for (@best) {
		    my %hit = %$_;
		    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
		    &update_pair($pair, \%hit, $aln, $i, $j, $term);
		}
		#die "Not implemented yet: unique for summative\n";
	    }
	}
    }

}

if ($mode eq "local") {
    for my $i (sort keys %$location) {
	my @sort = sort{
	    if ($a->{'from'} == $b->{'from'}) {
		$a->{'to'} <=> $b->{'to'};
	    } else {
		$a->{'from'} <=> $b->{'from'};
	    }
	} @{ $location->{$i} };
	my $a = shift @sort;
	# Create an array for equaly good hits that overlap
	my @equals;
	# Set tolerance ('e') for overlap
	my $e = $maxoverlap;
	# Update tolerance if percent option is selected
	$e = ($a->{'to'} - $a->{'from'} + 1) * $percentoverlap if $percentoverlap;
	# Select best (local)
	while (@sort) {
	    my $b = shift @sort;
	    # Check for overlap
	    #  Thanks to the sorting we need to check only that it starts before the other ends
	    my $tolerance = $e;
	    # Update tolerance if percent option is selected
	    ($tolerance) = &min( $e, ($b->{'to'} - $b->{'from'} + 1) * $percentoverlap ) if $percentoverlap;
	    if ($b->{'from'} + $tolerance <= $a->{'to'}) {
		if ($b->{'score'} > $a->{'score'}) {
		    # Keep B for now and "remove" A (and its equals)
		    $a = $b;
		    # Review equals and keep non-overlaping ones
		    my @keep;
		    for (@equals) {
			my $check = $tolerance;
			($check) = &min( $check, ($_->{'to'} - $_->{'from'} + 1) * $percentoverlap ) if $percentoverlap;
			unless ($b->{'from'} + $check <= $_->{'to'}) {
			    push @keep, $_;
			}
		    }
		    @equals = @keep;
		} elsif ($b->{'score'} == $a->{'score'}) {
		    # Store equaly good hits until N is reached if N is set
		    if ( ! $n || 1 + scalar(@equals) < $n ) {
			#unless ($n && 1 + scalar(@equals) <= $n) {
			push @equals, $a;
			$a = $b;
			# Update tolerance if percent option is selected
			$e = ($a->{'to'} - $a->{'from'} + 1) * $percentoverlap if $percentoverlap;
		    }
		    # Otherwise we stick with A and not use B
		}
		# Otherwise we stick with A and not use B
	    } else {
		# No overlap with A, so save A
		# Print equals and empty @equals
		while (@equals) {
		    my $x = shift @equals;
		    # print biointsam::sam_string($x->{'hit'}), "\n";
		    &print_sam($x->{'hit'});
		}
		# print biointsam::sam_string($a->{'hit'}), "\n";
		&print_sam($a->{'hit'});
		$a = $b;
		# Update tolerance if percent option is selected
		$e = ($a->{'to'} - $a->{'from'} + 1) * $percentoverlap if $percentoverlap;
	    }
	}
	while (@equals) {
	    my $x = shift @equals;
	    # print biointsam::sam_string($x->{'hit'}), "\n";
	    &print_sam($x->{'hit'});
	}
	# print biointsam::sam_string($a->{'hit'}), "\n" if $a;
	&print_sam($a->{'hit'});
    }    
} elsif ($mode eq "summative") {

    #print "# ", join("\t", "A-ID", "B-ID", "Similarity", "Identical", "Aligned", "A-coverage", "B-coverage", "A-covered", "B-covered", "A-length", "B-length", "orientation", "A-region", "B-region"), "\n";
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
	    $pair->{"PI:f"} = $pair->{"identical"} / $pair->{"aligned"};
	    $pair->{"max"} = $pair->{"identical"} / (1 + $pair->{"diff"});
	}
    }	
    
    my %found;
    
    # Best for each query
    my @best;
    for my $j (sort keys %query) {
	next unless $q2r->{$j};
	#    my $term = "PI:f"; # Used to be "similarity";
	#    my $term = "max"; # Not used currently
	#    my $term = "AS:i";
	#    my $term = "BS:i";
	# Select best (global)
	# Collect scores
	my %scores;
	#             sort based on specific score                                       register score in scores
	my @sorted = sort{ if ($q2r->{$j}->{$b}->{$term} == $q2r->{$j}->{$a}->{$term}) {
	                       $ref{$b} <=> $ref{$a}; # Prefere the longer reference when tied
			   } else {
			       $q2r->{$j}->{$b}->{$term} <=> $q2r->{$j}->{$a}->{$term};
			   }
                                                                                        } map { $scores{ $q2r->{$j}->{$_}->{$term} }++; $_} keys %{ $q2r->{$j} };
	#my @sorted = sort{$q2r->{$j}->{$b}->{$term} <=> $q2r->{$j}->{$a}->{$term}} map { $scores{ $q2r->{$j}->{$_}->{$term} }++; $_} keys %{ $q2r->{$j} };
	# sort scores to select top cutoff
	my @best_scores = sort{ $b <=> $a } keys %scores;
	# Get lowest score in the top N ($top)
	my $cutoff = $best_scores[$top - 1];
	my @best_set = grep{ $q2r->{$j}->{$_}->{$term} >= $cutoff } @sorted;
	
	#my ($i) = @sorted;
	my $count;
	for my $i (@best_set) {
	    $count++;
	    # Skip the rest if N is set and already reached
	    last if $n && $count > $n;
	    # Store info for the best one(s)
	    push @best, $q2r->{$j}->{$i};
	    $q2r->{$j}->{$i}->{"A"}->{"ID"} = $i;
	    $q2r->{$j}->{$i}->{"B"}->{"ID"} = $j;
	}
    }
    for my $pair (@best) {
	for my $hit (@{ $pair->{"hits"} }) {
	    # print biointsam::sam_string($hit), "\n";
	    &print_sam($hit);
	}
    }
}


#===SUBROUTINES=================================================================

sub print_sam {
    my ($hash) = @_;

    # Remove read pair flags ('/1' and '/2' for forward and reverse reads of a pair)
    if ($hash->{'FLAG'} & 64) {
	$hash->{'QNAME'} =~ s/\/1$//;
    } elsif ($hash->{'FLAG'} & 128) { 
	$hash->{'QNAME'} =~ s/\/2$//;
    }
    print biointsam::sam_string($hash), "\n";
}

sub update_pair {
    # Update pair info
    # This is used for summative mode
    my ($pair, $hitref, $aln, $i, $j, $term) = @_;
    my %hit = %$hitref;

    if ($term eq 'AS:i') {
	unless (defined $hit{'AS:i'}) {
	    # Make SAM output invalid to make the mistake obvious and die with error
	    print "ERROR: Aborted, due to missing alignment score ('AS:i:<int>' column). Check STDERR for more info.\n";
	    die "ERROR: Hit does not have an alignment score ('AS:i:<int>' column).\n" .
		"\tYou can use a tool to calculate alignment score before rerunning.\n" .
		"\tThe error occured in the following line:\n$_\n";
	}
	$pair->{$term} += $hit{$term};
    } elsif ($term eq 'BS:f') {
	unless (defined $hit{'BS:f'}) {
	    # Make SAM output invalid to make the mistake obvious and die with error
	    print "ERROR: Aborted, due to missing bitscore ('BS:f:<num>' column). Check STDERR for more info.\n";
	    die "ERROR: Hit does not have a bitscore ('BS:f:<num>' column).\n" .
		"\tMake sure to use the bitscore option only if you have a SAM file based on a BLAST run\n" .
		"\tThe error occured in the following line:\n$_\n";
	}
	$pair->{$term} += $hit{$term};
    }	

    my $sim = ($aln->{'length'} - $hit{'NM:i'}) / $aln->{'length'};
    $pair->{"identical"} += $aln->{'length'} - $hit{'NM:i'};
    $pair->{"diff"} += $hit{'NM:i'};
    $pair->{"aligned"} += $aln->{'length'};
    push @{ $pair->{"hits"} }, \%hit;
    push @{ $pair->{"similarity scores"} }, $sim;

    if ($term eq "PI:f") {
	$pair->{$term} = $pair->{"identical"} / $pair->{"aligned"} * 100;
    }
    
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
    # Done
}

sub print_help {
    # Print out the usage to STDERR
    # Takes in the @ARGV as input
    my ($args, $error) = @_;
    $error = "" unless $error;
    if ($error) {
	die $error . "\n$usage\n$description\n$options";
    }
    for (@$args) {
	if (/^-?-h(elp)?$/) {
	    print "$usage\n$description\n$options";
	    exit;
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

sub min {
    my (@array) = @_;
    # Return the element with the lowest numerical value
    my @sorted = sort{$a<=>$b} @array;
    return $sorted[0];
}
