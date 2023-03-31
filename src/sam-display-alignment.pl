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
    "A tool to display alignemnts of SAM files for exploratory purposes.\n" .
    "\tYou need to specify an input SAM file and the FASTA file for the reference.\n";
my $usage = 
    "Usage:\n\t$0 [-h] <SAM file> <Reference FASTA file>\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-w=<int> | --width=<int>\n\t\tSpecify the width of the coverage plots to be displayed (default is $covw).\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================


# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my @keep;
for (@ARGV) {
    if (/^--?w(?:idth)?=(\d+)$/) {
	$covw = $1;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;

my ($sam, $ref) = @ARGV;

# Print help if needed
biointbasics::print_help(\@ARGV, $info, 
			 "ERROR: you need to pass two files as arguments. " . 
			 "First the SAM file, then the reference FASTA file\n") unless $ref;



#print '%', ("123456789^" x 10), "\n";
my %ref;
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
    biointsam::parse_sam($_, \%ref, \%hit);
    # Header lines contain no hits
    unless (%hit) {
	# Print header to maintain a valid SAM output
#	print "$_\n";
	next;
    }
    # Nothing to do if there is no hit for the query sequence, so skip it
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";

    $hit{'SEQ'} = uc( $hit{'SEQ'} );
    
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
	    $hit{'SEQ'} = $seq;
	}
    } else {
	unless ($hit{'CIGAR'} =~ /H/) {
	    # Hard clipped sequences can generate errors, so skip them
	    $qname = $hit{'QNAME'};
	    $mem = $hit{'SEQ'};
	    $memflag = $hit{'FLAG'};
	}
    }

    print biointsam::sam_string(\%hit), "\n";

    print "\n"; 
    print "Coverage of the reference by the hit:\n";
    &plothit(\%hit, \%ref, $covw);
    print "\n";
    &displayhit(\%hit, $ref, $seqw);

    print "Coverage of the query by the hit:\n";
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    $query{ $hit{'QNAME'} } = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};
    biointsam::flip(\%hit, \%ref);
    &plothit(\%hit, \%query, $covw, 'flipped');
    print "\n";
    print "-" x ($seqw + 25), "\n";

    $hits++;
}
print "No hits to plot in the SAM file\n" unless $hits;

sub displayhit{
    # Displays the alignment of the hit in a layout similar to exonerate
    my ($href, $ref, $width) = @_;
    my %hit = %$href;
    if (! $hit{'SEQ'} || $hit{'SEQ'} eq "*") {
	print "\tUnable to print alignment! SEQ is absent for the entry\n\n";
	return;
    }
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});

    $width = 100 unless $width;
    my $rev;
    $rev++ if $hit{'FLAG'} & 16;
    
    # Get the sequence of R
    my $rio = Bio::SeqIO->new(-file=>$ref,-format=>'fasta');
    my $rseq;
    while ($rseq = $rio->next_seq()) {
	last if $rseq->primary_id() eq $hit{'RNAME'};
    }

    # Get R info
    my $start = $hit{'POS'};
    my $aligned = $aln->{'length'} - $aln->{'insertion'};
    my $rlen = $rseq->length();
    my $seq1 = substr($rseq->seq(), $start - 1, $aligned);

    # Get Q info
    my $seq2 = $hit{'SEQ'};
    my $qlen = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};


    # Add gaps based on CIGAR
    my $aligned1 = "";
    my $aligned2 = "";
    my $cigar = $hit{'CIGAR'};
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
    my $id1 = $hit{'RNAME'};
    my $id2 = $hit{'QNAME'};
    # collect match info and position mapping
    my $match;
    my $same = 0;
    my %p2seq1;
    my %p2seq2;
    my $p1 = $hit{'POS'} - 1;
    my $p2 = $aln->{'start'} - 1;
    my $d = 1;
    # Modify values for reverse complement
    if ($rev) {
	$d = -1;
	$p2 = $aln->{'start'} + $aln->{'length'} - $aln->{'deletion'};
    }
    for my $p (1..$total) {
	my $base1 = substr($seq1, $p - 1, 1);
	my $base2 = substr($seq2, $p - 1, 1);
	if ($base1 ne "-") {
	    $p1++;
	}
	if ($base2 ne "-") {
	    $p2 += $d;
	}
	$p2seq1{$p} = $p1;
	$p2seq2{$p} = $p2;
	my $m = &iupac_match($base1, $base2);
	$match .= $m;
	unless ($m eq ' ') {
	    $same++;
	}
    }
    my $pos = 0;
    my ($short1) = split/\s/, $id1;
    my ($short2) = split/\s/, $id2;
#    print "Title: $short1 vs $short2\n";
#    print "A: $id1\n";
#    print "B: $id2\n";
    my $sim = sprintf("%.2f", 100 * $same / $total);
    my $sim1 = sprintf("%.2f", 100 * $same / $rlen);
    my $sim2 = sprintf("%.2f", 100 * $same / $qlen);
    print "Alignment: $total; Identical: $same; Similarity: $sim\%;",
    #" A: $rlen ($sim1\%); B: $qlen ($sim2\%)",
    "\n";
    print "\n";
    while ($seq1) {
	my $w = $width;
	if (length $seq1 < $w) {
	    $w = length $seq1;
	}
	my $top = substr($seq1, 0, $w, '');
	my $mid = substr($match, 0, $w, '');
	my $bottom = substr($seq2, 0, $w, '');
	my $len = length $top;
	my $from = $pos + 1;
	my $to = $pos + $len;
	printf "%8s %s %8s\n", $p2seq1{$from}, $top, $p2seq1{$to};
	printf "%8s %s %8s\n", $from, $mid, $to;
	printf "%8s %s %8s\n", $p2seq2{$from}, $bottom, $p2seq2{$to};
	print "\n";
	$pos += $len;
    }
}

sub plothit{
    # Prints an ASCII plot for the hit to STDOUT
    my ($href, $ref, $width, $flipped) = @_;
    my %hit = %$href;
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});

    # Get R info
    my $start = $hit{'POS'};
    my $aligned = $aln->{'length'} - $aln->{'insertion'};
    my $len = $ref->{ $hit{'RNAME'} };

    # Get Q info
    my $q = $aln->{'length'} - $aln->{'deletion'};
    my $qlen = $q + $aln->{'unmapped'};

    # Calculate postions
    my $a;
    my $b;
    $a = sprintf("%.0f", ($start -1) / $len * $width);
    #$a = sprintf("%.0f", $start / $len * $width);
    $b = sprintf("%.0f", $aligned / $len * $width);
    # Add a point even if too small
    $b = 1 unless $b;
    my $c = $width - ($a + $b);
#    print "# \$a = $a; \$b = $b; \$c = $c; w = $width; aligned = $aligned\n" if ($a < 0 || $b < 0 || $c < 0);
    if ($c < 0) {
	# If 45.5 is gap (a) and 54.5 is match (b) then it will be rounded
	#  to 46 and 55. And that is longer than the width.
	# Adjsut so it still fits and give preference to match over gap.
	$a--;
	$c++;
    }

    # Calculate identity
    my $identity = ($aln->{'length'} - $hit{'NM:i'}) / $aln->{'length'};
    
    # Plot
    my $match = ">";
    $match = "<" if $hit{'FLAG'} & 16;
    my $gap = "_";
#    print "# \$a = $a; \$b = $b; \$c = $c\n" if ($a < 0 || $b < 0 || $c < 0);
    print "|" . ($gap x $a) . ($match x $b) . ($gap x $c) . "|";
    my $qblock = $hit{'QNAME'} . " (". sprintf("%.2f", $q / $qlen * 100) ."\%; l=" . $qlen . ")";
    my $rblock = $hit{'RNAME'} . " (". sprintf("%.2f", $aligned / $len * 100) ."\%; l=" . $len . ")";
    if ($flipped) {
	print " R:$qblock Q:$rblock";
    } else {
	print " R:$rblock Q:$qblock";
    }
    print " I:" . sprintf("%.4f", $identity) ."\n";
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
