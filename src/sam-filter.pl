#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

my $programname = "sam-filer.pl";
my $version = "1.1";
my $cmd = join(" ", $programname, @ARGV);

#===DESCRIPTION=================================================================

my $description = 
    "Description:\n\t" .
    "A tool to filter SAM files based on user specifications.\n" .
    "\tThe tool either opens the file specified as input or reads from STDIN when no file is given.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS]  [SAM file]\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-minlen=<int> | --minlen=<int>\n\t\tKeep only query sequences that at least <int> long.\n" .
    "\t-minaln=<int> | --minaln=<int>\n\t\tKeep only query sequences that have at least <int> long alignment with the reference.\n" .
    "\t-minsim=<float> | --minsim=<float>\n\t\tKeep only query sequences that have at least <float> similarity score.\n" .
    "\t-qcov=<float> | --query-coverage=<float>\n\t\tKeep only query sequences that have at least <float> part covered by the hit.\n" .
    "\t-rcov=<float> | --reference-coverage=<float>\n\t\tKeep only query sequences that cover at least <float> part of the reference sequence.\n" .
    "\t-mincomplexity=<float> | --mincomplexity=<float>\n\t\tKeep only hits where the aligned sequence has at least <float> 7-mer complexity (observed k-mers/maximal possible k-mers).\n" .
    "\t-dust | --dust=<int>\n\t\tKeep only hits where the aligned sequence has at most <int> (7 if not specified) dust score. (Algorithm adapted from PRINSEQ-lite).\n" .
    "\t-entropy | --entropy=<int>\n\t\tKeep only hits where the aligned sequence has at least <int> (70 if not specified) entropy score. (Algorithm adapted from PRINSEQ-lite).\n" .
    "\t-v\n\t\tInvert selection: report sequence that fail the requirement.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

my %ref;
my %query;

my $minlen;
my $minaln;
my $minsim;
my $mincomplexity;
my $dust;
my $entropy;
my $invert;
my $unmapped;
my $qcov;
my $rcov;
my @keep;

for (@ARGV) {
    if (/^--?minlen=(\d+)$/) {
	$minlen = $1;
    } elsif (/^--?q(?:uery-)?cov(?:erage)?=(\d\.\d+)$/) {
	$qcov = $1;
    } elsif (/^--?r(?:reference-)?cov(?:erage)?=(\d\.\d+)$/) {
	$rcov = $1;
    } elsif (/^--?minaln=(\d+)$/) {
	$minaln = $1;
    } elsif (/^--?mincomplexity=(\d\.\d+)$/) {
	$mincomplexity = $1;
    } elsif (/^--?dust(=\d+)?$/) {
	$dust = 7;
	$dust = $1 if /^--?dust=(\d+)$/;
    } elsif (/^--?entropy(=\d+)?$/) {
	$entropy = 7;
	$entropy = $1 if /^--?entropy=(\d+)$/;
    } elsif (/^--?minsim=(\S+)$/) {
	$minsim = $1;
    } elsif (/^-v$/) {
	$invert = "yes";
#	print STDERR "Using inverted mode: only those lines are reported that do not pass\n";
    } elsif (/^--?k(eep)?$/) {
	$unmapped = "yes";
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;

# Print help if needed
biointbasics::print_help(\@ARGV, $info);

# Remember read sequence until QNAME changes
my $qname = "";
my $mem = "";
my $memflag = 0;
my $header = "true";
my $previousprogram = "";
while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    unless (%hit) {
	print "$_\n";
	if (/^\@PG\tID:(\S+)/) {
	    $previousprogram = $1;
	}
	next;
    }
    if ($header) { # First line after the header section
	my $text = "\@PG\tID:$programname\tPN:$programname";
	$text .= "\tPP:$previousprogram" if $previousprogram;
	print $text . "\tVN:$version\tCL:$cmd\n";
	$header = undef;
    }
    
    if ($hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*") {
	print biointsam::sam_string(\%hit), "\n" if $unmapped;
	next;
    }

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
    
    # Get alginment data
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});

    # Apply filters if set
    my $failed;
    if ($minlen) {
	$failed++ if $minlen > $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};
    }
    if ($minaln && ! $failed) {
	$failed++ if $minaln > $aln->{'length'};
    }
    if ($minsim && ! $failed) {
	$failed++ if $minsim > ($aln->{'length'} - $hit{'NM:i'}) / $aln->{'length'};
    }
    if ($mincomplexity && ! $failed) {
	my @complexity = &kmer_complexity($aln->{'seq'});
	if ($complexity[0] eq "NA") {
	    print STDERR "WARNING: could not check k-mer complexity for $_\n";
	} else {
	    $hit{"kc:f"} = $complexity[0];
	    $failed++ if $mincomplexity > $complexity[0];
	}
    }
    if ($dust && ! $failed) {
	my $score = &dust_score($aln->{'seq'});
	if ($score eq "NA") {
	    print STDERR "WARNING: could not check DUST score for $_\n";
	} else {
	    $hit{"ds:i"} = $score;
	    $failed++ if $dust < $score;
	}
    }
    if ($entropy && ! $failed) {
	my $score = &entropy_score($aln->{'seq'});
	if ($score eq "NA") {
	    print STDERR "WARNING: could not check DUST score for $_\n";
	} else {
	    $hit{"ce:i"} = $score;
	    $failed++ if $entropy > $score;
	}
    }
    if ($qcov && ! $failed) {
	my $len = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};
	my $covered = $aln->{'length'} - $aln->{'deletion'};
	$failed++ if ($covered / $len) < $qcov;
    } 
    if ($rcov && ! $failed) {
	my $len = $ref{$hit{'RNAME'}};
	my $covered = $aln->{'length'} - $aln->{'insertion'};
	$failed++ if ($covered / $len) < $rcov;
    } 
    
    if ((! $failed && ! $invert) || ($failed && $invert)) {
	# Print if it passed and not inverted mode
	# Print if it did not pass but invert is selected
	print biointsam::sam_string(\%hit), "\n";
    }
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

sub dust_score {
    my ($seq) = @_;
    unless ($seq) {
	return "NA";
    }

    my $length = length($seq);
    
    my $WINDOWSIZE = 64;
    my $WINDOWSTEP = 32;
    my $WORDSIZE = 3;
    my @WINDOWSIZEARRAY = (0..61);
    my $POINTFIVE = 1/2;

    my $dustscore;
    my ($rest,$steps,@vals,$str,$num,$bynum);
    if($length <= $WINDOWSIZE) {
	$rest = $length;
	$steps = 0;
    } else {
	$steps = int(($length - $WINDOWSIZE) / $WINDOWSTEP) + 1;
	$rest = $length - $steps * $WINDOWSTEP;
	unless($rest > $WINDOWSTEP) {
	    $rest += $WINDOWSTEP;
	    $steps--;
	}
    }
    $num = $WINDOWSIZE-2;
    $bynum = 1/$num;
    $num--;
    my $mean = 0;
    foreach my $i (0..$steps-1) {
	$str = substr($seq,($i * $WINDOWSTEP),$WINDOWSIZE);
	my %counts = ();
	foreach my $i (@WINDOWSIZEARRAY) {
	    $counts{substr($str,$i,3)}++;
	}
	$dustscore = 0;
	foreach(values %counts) {
	    $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
	}
	push(@vals,($dustscore * $bynum));
    }
    #last step
    if($rest > 5) {
	$str = substr($seq,($steps * $WINDOWSTEP),$rest);
	my %counts = ();
	$num = $rest-2;
	foreach my $i (0..($num - 1)) {
	    $counts{substr($str,$i,3)}++;
	}
	$dustscore = 0;
	foreach(values %counts) {
	    $dustscore += ($_ * ($_ - 1) * $POINTFIVE);
	}
	push(@vals,(($dustscore / ($num-1)) * (($WINDOWSIZE - 2) / $num)));
    } else {
	push(@vals,31); #to assign a maximum score based on the scaling factor 100/31
    }
    $mean = &getArrayMean(@vals);
    my $score = int($mean * 100 / 31);
    return $score;

    sub getArrayMean {
	return @_ ? sum(@_) / @_ : 0;
    }
    sub sum {
	my $sum = 0;
	for (@_) {
	    $sum += $_;
	}
	return $sum;
    }
}

sub entropy_score {
    my ($seq) = @_;
    unless ($seq) {
	return "NA";
    }

    my $length = length($seq);
    
    my $WINDOWSIZE = 64;
    my $WINDOWSTEP = 32;
    my $WORDSIZE = 3;
    my @WINDOWSIZEARRAY = (0..61);
    my $LOG62 = log(62);
    my $ONEOVERLOG62 = 1/log(62);

    my ($rest,$steps,@vals,$str,$num,$bynum);
    if($length <= $WINDOWSIZE) {
	$rest = $length;
	$steps = 0;
    } else {
	$steps = int(($length - $WINDOWSIZE) / $WINDOWSTEP) + 1;
	$rest = $length - $steps * $WINDOWSTEP;
	unless($rest > $WINDOWSTEP) {
	    $rest += $WINDOWSTEP;
	    $steps--;
	}
    }
    $num = $WINDOWSIZE-2;
    $bynum = 1/$num;
    $num--;
    my $mean = 0;
    
    my $entropyval;
    foreach my $i (0..$steps-1) {
	$str = substr($seq,($i * $WINDOWSTEP),$WINDOWSIZE);
	my %counts = ();
	foreach my $i (@WINDOWSIZEARRAY) {
	    $counts{substr($str,$i,3)}++;
	}
	$entropyval = 0;
	foreach(values %counts) {
	    $entropyval -= ($_ * $bynum) * log($_ * $bynum);
	}
	push(@vals,($entropyval * $ONEOVERLOG62));
    }
    #last step
    if($rest > 5) {
	$str = substr($seq,($steps * $WINDOWSTEP),$rest);
	my %counts = ();
	$num = $rest-2;
	foreach my $i (0..($num - 1)) {
	    $counts{substr($str,$i,3)}++;
	}
	$entropyval = 0;
	$bynum = 1/$num;
	foreach(values %counts) {
	    $entropyval -= ($_ * $bynum) * log($_ * $bynum);
	}
	push(@vals,($entropyval / log($num)));
    } else {
	push(@vals,0);
    }
    $mean = &getArrayMean(@vals);
    return int($mean * 100);
}
