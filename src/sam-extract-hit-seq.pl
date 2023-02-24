#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

use Bio::SeqIO;
use Bio::Seq;

my $covw = 80;
my $seqw = 60;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "This script extracts hit (target, query, read) seqequences from SAM files\n" .
    "\twith some extra information:\n" .
    "\t - location in query sequence\n" .
    "\t - similarity (identical positions / alignment length)\n" .
    "\t - coverage breadth of the reference sequence\n" .
    "\t - unmapped positions of the refernce sequence at the two ends\n";
my $usage = 
    "Usage:\n\t$0 [-h] <SAM file> <Reference FASTA file>\n";
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

my %ref;
while(<>) {
    my %hit;
    &parse_sam($_, \%ref, \%hit); 
    next unless %hit;
    next if $hit{'RNAME'} eq "*";
    # Get alginment data
    my $aln = &parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    my $len = $aln->{'length'};
    my $del = $aln->{'deletion'};
    my $ins = $aln->{'insertion'};
    my $sim = sprintf "%.2f", ( $len - $hit{'NM:i'} ) / $len;
    my $cov = sprintf "%.2f", ( $len - $ins ) / $ref{ $hit{'RNAME'} };
    # Is something missing from reference
    my $mis5 = $hit{'POS'} - 1;
    my $mis3 = $ref{ $hit{'RNAME'} } - ($len - $ins) - $mis5;
    # Get query region position
    my $pos;
    my $qlen = $len - $del;
    $pos = $aln->{'start'} . ".." . ($aln->{'start'} - 1 + $qlen);
    $pos = "complement(" . $pos . ")" if $aln->{'reverse'};
    my $id = $hit{'QNAME'} . ":" . $pos;
    my $sep = "\t";
    my $def = join($sep, "ref:" . $hit{'RNAME'}, "sim:$sim", "cov:$cov", "len:$qlen", "missing:($mis5;$mis3)");
    if ($hit{'SEQ'} eq '*') {
	print STDERR "Warning: skipping line ($hit{'QNAME'}), because seq is '*'\n";
	next;
    }
    print &to_fasta($id . $sep . $def, $aln->{'seq'});
}

sub parse_sam {
    # Reads a SAM file line by line and converts each entry to a hash
    my ($line, $ref, $hash) = @_;
    my @field = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;    
    s/\R+//;
    my @col = split/\t/;
    if (/\@/) { # header
	my %h;
	for (@col[1..$#col]) {
	    /^([^:]+):(.*)/;
	    $h{$1} = $2;
	}
	if ($col[0] eq '@SQ') {
	    $ref->{ $h{'SN'} } = $h{'LN'};
	}
    } else {
	my $i = 0;
	for (@col) {
	    if ($i < scalar @field) {
		$hash->{ $field[$i] } = $col[$i];
	    } else {
		/^([^:]+\:[^:]+):(.*)/;
		$hash->{$1} = $2;
	    }
	    $i++;
	}
    }
    return;
}

sub to_fasta {
    # Return a fasta formated string
    my ($seq_name, $seq, $len) = @_;
    # default to 60 characters of sequence per line
    $len = 60 unless $len;
    # Print ID line
    my $formatted_seq = ">$seq_name\n";
    # Add sequence lines with $len length
    while (my $chunk = substr($seq, 0, $len, "")) {
	$formatted_seq .= "$chunk\n";
    }
    return $formatted_seq;
}

sub parse_cigar {
    my ($cigar, $bitflag, $seq) = @_;
    $seq = undef if $seq eq '*';
    my $len = 0;
    my $insert = 0;
    my $deletion = 0;
    my $start = 1;
    my $unmapped = 0;
    # Remove clipped ends and store info
    if ($cigar =~ s/^(\d+)([SH])//) {
	my $length = $1;
	$start += $length unless $bitflag & 16;
	$unmapped += $length;
	if ($2 eq 'S' && $seq) {
	    $seq =~ s/^.{$length}//;
	}
    }
    if ($cigar =~ s/(\d+)([SH])$//) {
	my $length = $1;
	$start += $length if $bitflag & 16;
	$unmapped += $length;
	if ($2 eq 'S' && $seq) {
	    $seq =~ s/.{$length}$//;
	}
    }
    # Parse element by element
    while (length $cigar) {
	$cigar =~ s/^(\d+)//;
	my $length = $1;
	$cigar =~ s/^(\D)//;
	my $type = $1;
	if ($type eq "D" || $type eq "N") {
	    $deletion += $length;
	    $len += $length;
	} elsif ($type eq "I") {
	    $len += $length;
	    $insert += $length;
	} elsif ($type eq "M" || $type eq "=" || $type eq "X") {
	    $len += $length;
	} elsif ($type eq "S" || $type eq "H") {
	    die "ERROR: unexpected $type element in cigar\n";
	} else {
	    print STDERR "Unexpected cigar element found ('$length$type') that will be skipped\n";
	}
    }
    return {"start" => $start,
	    "length" => $len,
	    "insertion" => $insert,
	    "deletion" => $deletion,
	    "unmapped" => $unmapped,
	    "reverse" => $bitflag & 16, # 0 or 16
	    "seq" => $seq,
    };
}
