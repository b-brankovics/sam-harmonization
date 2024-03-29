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
    "A tool to calculate ANI and coverage statistics for SAM files.\n" .
    "\tThe tool either opens the file specified as input or reads from STDIN when no file is given.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] [SAM file] [Reference FASTA] [Query FASTA]\n";
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

# Get input files
# - none => STDIN
# - one => only SAM file
# - three
my ($sam, $ref_fasta, $query_fasta) = @ARGV;

# Check input files
my $error = "";
if (@ARGV) {
    # At least SAM file is specified
    unless($sam ne '-' && -e $sam && ! -z $sam){
            $error .= "ERROR: SAM file ('$sam') is empty or missing\n";
    }  
    if (@ARGV > 3) {
        $error .= "ERROR: Too many inputs specified, please use maximum 3.\n";
    } elsif (@ARGV == 2) {
        $error .= "ERROR: If you wish to specify FASTA files, please specify both reference and query files.\n";
    } elsif (@ARGV == 3) {
        # Check the FASTA files
        if (! -e $ref_fasta || -z $ref_fasta) {
            $error .= "ERROR: Reference FASTA file ('$ref_fasta') is empty or missing\n";
        }
        if (! -e $query_fasta || -z $query_fasta) {
            $error .= "ERROR: Query FASTA file ('$query_fasta') is empty or missing\n";
        }
    }
}
# Print help or errors if needed
biointbasics::print_help(\@ARGV, $info, $error);

# Store length info based on SAM file for IDs
my %ref;
my %query;

# ANI calculation
my $numer; # nummerator
my $denom; # denominator
my $length;
my $mis;
# Store homology stretches per input to use for coverage calculation
my %r;
my %q;
# Keep all the individual identity scores
my @scores;

my $samfh;
# Read from STDIN if no file is specified or '-' is used
if (! $sam || $sam eq '-') {
    $samfh = *STDIN;
} else {
    open($samfh, '<', $sam) || die $!;
}
while(<$samfh>) {
    # Parse SAM line by line and populate %ref from the header and save hit info
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    # Next unless there is hit or mapping
    next unless %hit;
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
    # Get alginment data for hit
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    # Store the length of query seq
    $query{ $hit{'QNAME'} } = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};
    # Alignment length and mismatches
    my $length = $aln->{'length'};
    my $mis = $hit{'NM:i'};

    # Update values for ANI calculation
    $numer += $length - $mis;
    $denom += $length;
    # Save individual score value
    push @scores, ($length - $mis) / $length;

    # store info for coverage
    # ref
    # Get start and end positions
    my $start1 = $hit{'POS'};
    my $end1 = $hit{'POS'} + $length - $aln->{'insertion'} - 1;
    &add_hit(\%r, $hit{'RNAME'}, $start1, $end1);
    # query
    # Get start and end positions
    my $start2 = $aln->{'start'};
    my $end2 = $aln->{'start'} - 1 + $aln->{'length'} - $aln->{'deletion'};
    &add_hit(\%q, $hit{'QNAME'}, $start2, $end2);
}

# Merge overlapping or continuous hits before calculating (breadth) coverage
&merge_overlapping(\%r);
&merge_overlapping(\%q);
my $ref_cov = &calc_coverage(\%r);
my $query_cov = &calc_coverage(\%q);

# Report ANI if there were any hits => denominator > 0
unless ($denom) {
    print "No alignment in the SAM file\n";
    exit;
}
print "ANI: " . (sprintf '%.2f', 100 * $numer / $denom) . "% ($numer/$denom)\n";

# Print score range and mean for indivudal alignments/hits
my @sort = sort {$a<=>$b} @scores;
my $sum;
for (@sort) {
    $sum += $_;
}
print "(Range of similarity scores [$sort[0], $sort[-1]]; arithmetic mean " . ($sum / scalar @sort) . ")\n";

# Get full lengths for coverage info
my ($ref_len, $query_len);
# Use FASTA files if possible
if ($ref_fasta && $query_fasta) {
    # store seq info hash: ID -> sequence
    my $ref_seq = {};
    my $query_seq = {};
    biointbasics::read_fasta($ref_seq, [], $ref_fasta);
    biointbasics::read_fasta($query_seq, [], $query_fasta);
    # Add up sequence lengths
    for (values %$ref_seq) {
        $ref_len += length($_);
    }
    for (values %$query_seq) {
        $query_len += length($_);
    }
} else {
    # Get length info from SAM data
    for (values %ref) {
        $ref_len += $_;
    }
    for (values %query) {
        $query_len += $_;
    }
    # Warn user
    print STDERR "WARNING: SAM file is used to calculate reference and query length. This may be inacurate if it was filtered.\n";
}

# Print coverage info
print "R covered: $ref_cov";
print " (", (sprintf '%.2f', (100 * $ref_cov / $ref_len)), "%)";
print "\n";
print "Q covered: $query_cov";
print " (", (sprintf '%.2f', (100 * $query_cov / $query_len)), "%)";
print "\n";

# Print summary in 2 line TSV
print "# ", join("\t", "Similarity", "Identical", "Aligned", "A-covered", "B-covered", "A-length", "B-length"), "\n";
print join("\t", $numer / $denom, $numer, $denom, $ref_cov, $query_cov, $ref_len, $query_len), "\n";

#===SUBROUTINES=================================================================

sub add_hit {
    # Add hit info to hash based on seq_ID, start and end position
    my ($hashref, $id, $start, $end) = @_;
    if ($hashref->{$id} && $hashref->{$id}->{$start}) {
        # if already exists a hit with the same start, extend if longer
	    if ($hashref->{$id}->{$start} < $end) {
	        $hashref->{$id}->{$start} = $end;
	    }
    } else {
        # Create new hit
	    $hashref->{$id}->{$start} = $end;
    }
}

sub merge_overlapping {
    # Merge overlaps or continuous stretches for ref or query hash
    # Modifies the hash content
    my ($hasref) = @_;
    for my $id (sort keys %$hasref) {
        my $from;
        my $to;
        for my $s (sort{$a <=> $b} keys %{ $hasref->{$id} }) {
    	    if ($to && $to + 1 >= $s) {
	            if ($hasref->{$id}->{$s} > $to) {
		        # extend the old one
		        $hasref->{$id}->{$from} = $hasref->{$id}->{$s};
	            }
	            # remove the current one
	            delete $hasref->{$id}->{$s};
	            $s = $from;
	        }
	        $from = $s;
	        $to = $hasref->{$id}->{$s};
        }
    }
}

sub calc_coverage {
    # Sum total coverage for ref or query hash and return coverage value
    my ($hasref) = @_;
    my $cov;
    for my $id (sort keys %$hasref) {
        for my $s (sort{$a <=> $b} keys %{ $hasref->{$id} }) {
	    $cov += $hasref->{$id}->{$s} + 1 - $s;
        }
    }
    return $cov;
}