#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

use Bio::SeqIO;
use Bio::Seq;

use Data::Dump qw(dump);

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to update the SEQ portion of a SAM file based on a FASTA or FASTQ file.\n" .
    "\tYou need to specify an input SAM file and the FAST[AQ] file for the query.\n";
my $usage = 
    "Usage:\n\t$0 [-h] <SAM file> <Query FASTQ file>\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-fq | --fastq\n\t\tSpecify that the query file is in FASTQ format (default is FASTA).\n" .
    "\t-r=<Reference FASTA file> | --ref=<Reference FASTA file>\n\t\tSpecify the reference FASTA file to update the header, and the alignment data if '--cigar' is selected.\n" .
    "\t-c | --cigar\n\t\tUpdate alignment data: CIGAR and NM:i:<int> (for the alignment score use 'sam-score.pl'.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================


# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my $rfile;
my $update;
my $fastq;

my @keep;
for (@ARGV) {
    if (/^--?r(?:ef)?=(\S+)$/) {
	$rfile = $1;
    } elsif (/^--?f(ast)?q$/) {
	$fastq++;
    } elsif (/^--?c(igar)?$/) {
	$update++;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;

my ($sam, $query) = @ARGV;

# Print help if needed
biointbasics::print_help(\@ARGV, $info, 
			 "ERROR: you need to pass two files as arguments. " . 
			 "First the SAM file, then the query FASTQ file\n") unless $query;

# Open Query file for reading
my $check;
if ($query eq '-') {
    $check = *STDIN;
} else {
    open ($check, '<', $query) || die $!;
}
my $firstline = <$check>;
if ($firstline =~ /^>/) {
    print STDERR "Warning: FASTQ was invoked, but file is FASTA format. Switching to FASTA format.\n" if $fastq;
    $fastq = undef;
} elsif ($firstline =~ /^\@/){
    print STDERR "Warning: FASTQ was not invoked, but file is FASTQ format. Switching to FASTQ format.\n" unless $fastq;
    $fastq++;
} else {
    print STDERR "Warning: Query input format verification failed.\n" unless $fastq;
}

my $qio;
if ($fastq) {
    $qio = Bio::SeqIO->new(-file=>$query,-format=>'fastq');
} else {
    $qio = Bio::SeqIO->new(-file=>$query,-format=>'fasta');
}




#print '%', ("123456789^" x 10), "\n";
my %ref;
# Update REF header
if ($rfile) {
    my $rio = Bio::SeqIO->new(-file=>$rfile,-format=>'fasta');
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
my $qseq;
while(<$samfh>) {
    #print "$_";
    my %hit;
    next if /^\@SQ/ && $rfile;
    biointsam::parse_sam($_, \%ref, \%hit);
    # Header lines contain no hits
    unless (%hit) {
	# Print header to maintain a valid SAM output
	print "$_\n";
	next;
    }
    unless ($qname eq $hit{'QNAME'}) {
	# Look up seq in FASTQ
	while ($qseq = $qio->next_seq()) {
	    last if $qseq->id() =~ /^$hit{'QNAME'}\b/;
	}
	# Save QNAME for the next entry
	$qname = $hit{'QNAME'};
    }

    my $start = 1;
    my $qlen;
    my $seq = $qseq->seq();
    my $qual = "*";
    if ($fastq) {
	# Write the whole FASTQ entry to the string
	open my $fh, ">", \$qual;
	my $qualOut=Bio::SeqIO->new(-fh=>$fh,-format=>"fastq");
	$qualOut->write_seq($qseq);
	# We need only the Quality string (last line), and no line endings
	$qual =~ /\+\R+(.+)/;
	$qual = $1;
    }

    unless($hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*") {

	&updatehit(\%hit, $rfile) if $update;
	
	if ($hit{'FLAG'} & 16) {
	    $seq = biointsam::reverse_seq($seq);
	    $qual = reverse $qual;
	}
	
	my $trash;
	if ($hit{'CIGAR'} =~ /^(\d+)H/) {
	    # clip sequence
	    $trash = substr($seq, 0, $1, '');
	    $trash = substr($qual, 0, $1, '');
	}
	if ($hit{'CIGAR'} =~ /(\d+)H$/) {
	    # clip sequence
	    my $l = $1;
	    $seq =~ s/.{$l}$//;
	    $qual =~ s/.{$l}$//;
	    #$trash = substr($seq, 0 - $1, '');
	    #$trash = substr($qual, 0 - $1, '');
	}
	#$hit{'SEQ'} = substr($seq, $start - 1, $qlen);
	$hit{'SEQ'} = $seq;

	$hit{'QUAL'} = $qual;
    } else {
	$hit{'SEQ'} = $seq;
	$hit{'QUAL'} = $qual;
    }
    
    print biointsam::sam_string(\%hit), "\n";
}


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

