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
    "A tool to print TSV format hit info for SAM files.\n" .
    "\tThe tool either opens the file specified as input or reads from STDIN when no file is given.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS]  [SAM file]\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "Output:\n" .
    "\tThe output is tab separate file with the following columns:\n" .
    "\t  1. Reference ID\n" .
    "\t  2. Query ID\n" .
    "\t  3. Sequence similarty (%)\n" .
    "\t  4. Hit start (both positions refer to the reference sequence)\n" .
    "\t  5. Hit end (if start > end, then it is reverse hit)\n" .
    "\t  6. 7-mer complexity of the aligned query sequence\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

my %ref;
my %query;
my $hits;

biointbasics::print_help(\@ARGV, $info);

# Remember read sequence until QNAME changes
my $qname = "";
my $mem = "";
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
    # Get seq in case it is a secondary alignment (only works if sorted by read name)
    if ($qname eq $hit{'QNAME'}) {
	$hit{'SEQ'} = $mem;
    } else {
	$qname = $hit{'QNAME'};
	$mem = $hit{'SEQ'};
    }

    # Nothing to do if there is no hit for the query sequence, so skip it
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
    &hit2tsv(\%hit);
    print "\n";
    $hits++;
}
print "No hits to plot in the SAM file\n" unless $hits;

sub hit2tsv{
    # Prints an ASCII plot for the hit to STDOUT
    my ($href) = @_;
    my %hit = %$href;
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});

    # Get R info
    my $start = $hit{'POS'};
    my $aligned = $aln->{'length'} - $aln->{'insertion'};
    my $stop = $start + $aligned - 1;

    ($start, $stop) = ($stop, $start) if $hit{'FLAG'} & 16;
    
    my $identity = ($aln->{'length'} - $hit{'NM:i'}) / $aln->{'length'};
    my @complexity = &kmer_complexity($aln->{'seq'});

    print join("\t", $hit{'RNAME'}, $hit{'QNAME'}, sprintf("%.4f", $identity) * 100, $start, $stop, $complexity[0]);
}


sub kmer_complexity {
    my ($seq, $k, $circular) = @_;

    $k = 7 unless $k;
    my $absolute_max = 4**$k;
    unless ($seq) {
	return ("NA", "NA", "NA");
    }
    my $max = length($seq) * 2;
    unless ($circular) {
	$max = (length($seq) - $k + 1) * 2;
    }
    my %kmers;
    &generate_kmers($seq, $k, \%kmers, $circular);
    my $count = scalar( keys %kmers );
    if ($max > $absolute_max) {
	$max = $absolute_max;
    }
    return ($count/$max, $count, $max);
}

sub generate_kmers {
    # generates kmers from a sequence
    my ($seq, $k, $kmer_ref, $circular) = @_;

    my $len = length $seq;
    # All kmers that fit
    for (0..($len-$k)) {
	my $kmer = substr $seq, $_, $k;
	$$kmer_ref{$kmer}++;
	$$kmer_ref{ &reverse_seq($kmer) }++;
    }
    if ($circular) {
	# Add circular kmers
	for (1..($k -1)) {
	    my $kmer = substr $seq, ($len - $k + $_);
	    $kmer .= substr $seq, 0, $_;
	    $$kmer_ref{$kmer}++;
	    $$kmer_ref{ &reverse_seq($kmer) }++;
	}
    }
    return;
}

sub reverse_seq {
    # Reverse complements the sequences
    my ($seq) = @_;
    # Reverse the sequnce
    my $complement = reverse $seq;
    # Complement the sequence
    $complement =~ tr/ACGTacgtWwMmRrSsKkYyBbVvDdHh/TGCAtgcaWwKkYySsMmRrVvBbHhDd/;
    return $complement;
}
