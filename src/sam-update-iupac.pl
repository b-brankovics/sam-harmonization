#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

use Bio::SeqIO;
use Bio::Seq;

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to update alignemnt info of SAM files for ambiguity sites in the reference.\n" .
    "\tYou need to specify an input SAM file and the FASTA file for the reference.\n";
my $usage = 
    "Usage:\n\t$0 [-h] <SAM file> <Reference FASTA file>\n";
my $options = 
    "Options:\n" .
    "\t-head | --header\n\t\tUpdate the reference sequence information in the header.\n" .
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

my $header;
my @keep;
for (@ARGV) {
    if (/^--?head(er)?$/) {
	$header++;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;

my ($sam, $ref) = @ARGV;

# Check that input exists
if ($ref) {
    if (! -e $ref || -z $ref) {
	biointbasics::print_help(\@ARGV, $info, "ERROR: '$ref' input reference FASTA does not exist or it is empty.\n");
    }
} else {
    biointbasics::print_help(\@ARGV, $info, 
			     "ERROR: you need to pass two files as arguments. " . 
			     "First the SAM file, then the reference FASTA file\n");
}



#print '%', ("123456789^" x 10), "\n";
my %ref;
# Update header if asked
if ($header) {
    my $rio = Bio::SeqIO->new(-file=>$ref,-format=>'fasta');
    while ( my $rseq = $rio->next_seq() ) {
	my $id = $rseq->id();
	my $len = $rseq->length();
	$ref{$id} = $len;
	print join("\t", '@SQ', 'SN:' . $id, 'LN:' . $len), "\n";
    }
}

my %query;
my $hits;
my $samfh;
if ($sam eq '-') {
    $samfh = *STDIN;
} else {
    open($samfh, '<', $sam) || die $!;
}
# Remember read sequence until QNAME changes
my $qname = "";
my $mem = "";
my $memflag = 0;
while(<$samfh>) {
    #print "$_";
    my %hit;
    next if /^\@SQ/ && $header;
    biointsam::parse_sam($_, \%ref, \%hit);
    # Header lines contain no hits
    unless (%hit) {
	# Print header to maintain a valid SAM output
	print "$_\n";
	next;
    }
    # Nothing to do if there is no hit for the query sequence, so skip it
    next if $hit{'RNAME'} eq "*";

    $hit{'SEQ'} = uc($hit{'SEQ'});
    
    # Get seq in case it is a secondary alignment (only works if sorted by read name)
    if ($qname eq $hit{'QNAME'}) {
	if ($hit{'SEQ'} eq '*') {
	    my $seq = $mem;
	    # reverse if orientation differs
	    $seq = biointsam::reverse_seq($seq) if ($hit{'FLAG'} & 16) != ($memflag & 16);
	    # Process hard clipping
	    if ($hit{'CIGAR'} =~ /^(\d+)H/) {
		$seq = substr($seq, $1);
	    }
	    if ($hit{'CIGAR'} =~ /(\d+)H$/) {
		$seq = reverse substr(reverse($seq), $1);
	    }
	    $hit{'SEQ'} = $mem;
	}
    } else {
	unless ($hit{'CIGAR'} =~ /H/) {
	    # Hard clipped sequences can generate errors, so skip them
	    $qname = $hit{'QNAME'};
	    $mem = $hit{'SEQ'};
	    $memflag = $hit{'FLAG'};
	}
    }

    &updatehit(\%hit, $ref);
    
    print biointsam::sam_string(\%hit), "\n";

    $hits++;
}

print "No hits to plot in the SAM file\n" unless $hits;

sub updatehit{
    # Update the alignment info for ambiguity sites
    my ($hit, $ref) = @_;
    if (! $hit->{'SEQ'} || $hit->{'SEQ'} eq "*") {
	print "\tUnable to print alignment! SEQ is absent for the entry\n\n";
	return;
    }
    my $aln = biointsam::parse_cigar($hit->{'CIGAR'}, $hit->{'FLAG'}, $hit->{'SEQ'});

    my $rev;
    $rev++ if $hit->{'FLAG'} & 16;
    
    # Get the sequence of R
    my $rio = Bio::SeqIO->new(-file=>$ref,-format=>'fasta');
    my $rseq;
    while ($rseq = $rio->next_seq()) {
	last if $rseq->primary_id() eq $hit->{'RNAME'};
    }

    # Get R info
    my $start = $hit->{'POS'};
    my $aligned = $aln->{'length'} - $aln->{'insertion'};
    my $rlen = $rseq->length();
    my $seq1 = substr($rseq->seq(), $start - 1, $aligned);

    # Get Q info
    my $seq2 = uc($hit->{'SEQ'});
    my $qlen = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};

    &update_cigar($hit);
 
    # Add gaps based on CIGAR
    my $aligned1 = "";
    my $aligned2 = "";
    my $cigar = $hit->{'CIGAR'};
    while ($cigar) {
	$cigar =~ s/^(\d+)(\D)//;
	my $len = $1;
	my $type = $2;
	# Skip hard clipping
	next if ($type eq "H");
	# Clip soft clipping
	if ($type eq "S") {
	    substr($seq2, 0, $len, "");
	} elsif ($type eq "M" || $type eq "=" || $type eq "X") {
	    $aligned1 .= substr($seq1, 0, $len, "");
	    $aligned2 .= substr($seq2, 0, $len, "");
	} elsif ($type eq "D" || $type eq "N") {
	    $aligned1 .= substr($seq1, 0, $len, "");
	    $aligned2 .= "-" x $len;
	} elsif ($type eq "I") {
	    $aligned1 .= "-" x $len;
	    $aligned2 .= substr($seq2, 0, $len, "");
	}
    }
    $seq1 = $aligned1;
    $seq2 = $aligned2;
    my $total = length $seq1;
    my $mismatches = 0;
    for my $p (1..$total) {
	my $base1 = substr($seq1, $p - 1, 1);
	my $base2 = substr($seq2, $p - 1, 1);
	if ($base1 eq "-" || $base2 eq "-") {
	    next;
	}
	my $m = &iupac_match($base1, $base2);
	if ($m eq ' ') {
	    $mismatches++;
	}
    }
    $hit->{'NM:i'} = $mismatches;
}

sub update_cigar {
    my ($hit) = @_;
    # Update CIGAR if needed
    my $len = 0;
    my $type;
    my $cigar = "";
    my $original = $hit->{'CIGAR'};
    while (length $original) {
	$original =~ s/^(\d+)//;
	my $l = $1;
	$original =~ s/^(\D)//;
	my $t = $1;
	if ($t eq '=' || $t eq 'X') {
	    $t = 'M';
	}
	if ($type) {
	    if ($t eq $type) {
		$len += $l;
	    } else {
		$cigar .= $len . $type;
		$len = $l;
		$type = $t;
	    }
	} else {
	    $len = $l;
	    $type = $t;
	}
    }
    $cigar .= $len . $type;
    $hit->{'CIGAR'} = $cigar;
}


sub iupac_match {
    # Check match to IUPAC code
    my ($r, $q) = @_;
    # IUPAC
    my %nt = (
	W => "[AT]",
	S => "[CG]",
	M => "[AC]",
	K => "[GT]",
	R => "[AG]",
	Y => "[CT]",
	B => "[CGT]",
	D => "[AGT]",
	H => "[ACT]",
	V => "[ACG]",
	N => "[ACGT]",
	);
    if ($nt{$r}) {
	if ($q =~ /$nt{$r}/) {
	    return ':';
	} else {
	    return ' ';
	}
    } else {
	if ($r eq $q) {
	    return '|';
	} else {
	    return ' ';
	}
    }
}
