#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointbasics;

#===DESCRIPTION=================================================================
# A tool to remove duplicate sequence from FASTA format files and
#  print the groups to STDERR

my $description = 
    "Description:\n\t" .
    "A tool to convert (nucmer) delta files to SAM format printed on STDOUT.\n" . 
    "\tSequence information needs to be read from the source FASTA files\n" .
    "\tspecified in the header of the delta file.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] [delta file]\n";
my $options = 
    "Options:\n" .
    "\t-h | --help\n\t\tPrint the help message; ignore other arguments.\n" .
    "\t-s | --skip\n\t\tSkip parsing the source FASTA files. The SAM file will contain\n" .
    "\t\tonly those sequence IDs that are present in the delta files\n" .
    "\t\t(have hits), and the SEQ column will contain '*' instead of the\n" .
    "\t\tactual sequence of the query.\n" .
    "\n";
my $info = {
    description => $description,
    usage => $usage,
    options => $options,
};

#===MAIN========================================================================

my $skip_files;
my @keep;
for (@ARGV) {
    if (/^--?s(kip)?$/) {
	$skip_files++;
    } else {
	push @keep, $_;
    }
}
@ARGV = @keep;


# Print help if needed
biointbasics::print_help(\@ARGV, $info);




my $mapq = 255;

my $ref;
my $query;

my %reflen;
my $r;
my $t;
my @mem;
my $cigar;
my $bitflag;
my $pos;
my $length;
my $covered;
my $tlen;
my $edit;
my $seq = "";
my $qual = "*";
my $head = "";
my $body = "";
my $start;
my $end;

my $fasta;
my %seen;
my $pg = "";
while(<>) {
    if ($. == 1) {
	# Store file names and their path
	s/\R+//;
	($ref, $query) = split/\s+/;
	unless ($skip_files) {
	    my $reffasta = &read_fasta($ref);
	    for (keys %$reffasta) {
		$reflen{$_} = length $reffasta->{$_};
	    }
	    $fasta = &read_fasta($query);
	}
	$pg = "\@PG\tID:nucmer\tPN:nucmer\tVN:3.1\tCL:nucmer $ref $query\n";
    }
    if (/^>(\S+) (\S+) (\d+) (\d+)/) {
	# Store sequence IDs
	$r = $1;
	$t = $2;
	$reflen{$r} = $3;
	$tlen = $4;
	$seen{$t}++;
    } elsif (/^(\d+) (\d+) (\d+) (\d+) (\d+)/) {
	$bitflag = 0;
	if ( ($1 > $2 && $3 < $4) || ($1 < $2 && $3 > $4) ) {
	    # Reverse match
	    $bitflag += 16;
	}
	#start1 end1  s2  e2  mismatches (including gaps)
	my ($s, $e) = sort {$a<=>$b} ($1, $2);
	my ($s2, $e2) = sort {$a<=>$b} ($3, $4);
	$edit = $5;
	$pos = $s;
	$length = $e2 - $s2 + 1;
	$covered = 0;
	# Process clipping
	my $clip1 = $s2 - 1;
	my $clip2 = $tlen - $e2;
	$start = "";
	$end = "";
	$start = $clip1 . "H" if $clip1;
	$end = $clip2 . "H" if $clip2;
	unless ($skip_files) {
	    $seq = substr($fasta->{$t}, $s2 - 1, $e2 - $s2 + 1);
	    $seq = &reverse_seq($seq) if $bitflag & 16;
	} else {
	    $seq = "*";
	}
	($start, $end) = ($end, $start) if $bitflag & 16;
    } elsif ($pos && /^(-?(\d+))/) {
	my $full = $1;
	my $l = $2;
	# Check for gaps in the alignemnt
	if ($full == 0) {
	    # End of the alignemnt, time to print CIGAR
	    if ($length > $covered) {
		push @mem, {'type' => 'M', 'length' => $length - $covered};
	    }
	    my $cigar = $start;
	    #my $cigar;
	    while (@mem) {
		$_ = shift @mem;
		$cigar .= $_->{'length'} . $_->{'type'};
	    }
	    # Add remaining match
	    
	    $cigar .= $end;
	    # Calculate alignment score (exonerate)
	    # Match     +5
	    # Mismatch  -4
	    # Gapextend -4
	    # Gapopen   -12 (all edit distance is -4 plus an extra -8 for each indel [DI])
	    my $gapcount = () = $cigar =~ /[DIN]/g;
	    my @matches = $cigar =~ /(\d+)M/g;
	    my @indel = $cigar =~ /(\d+)[DIN]/g;
	    my $sum = 0;
	    for (@matches) {
		$sum += $_;
	    }
	    my $gappos = 0;
	    for (@indel) {
		$gappos += $_;
	    }
	    my $identical = $sum - ($edit - $gappos);
	    my $score = 5 * $identical - 4 * $edit - 8 * $gapcount;
	    $body .= join("\t", $t, $bitflag, $r, $pos, $mapq, $cigar, "*", 0, 0, $seq, $qual, "NM:i:" . $edit, "AS:i:" . $score) . "\n";
	    
	    $pos = undef;
	} elsif ($full < 0) {
	    $covered += 1;
	    # gap added to the reference => I (already included in the $mis)
	    if (@mem && $mem[-1]->{'type'} eq "I" && $l == 1) {
		# increment previous
		$mem[-1]->{'length'} += 1;
	    } else {
		if ($l > 1) {
		    push @mem, {'type' => 'M', 'length' => $l - 1};
		    $covered += $l - 1; 
		}
		push @mem, {'type' => 'I', 'length' => 1};
	    }
	} elsif ($full > 0 ) {
	    # gaps in query => D
	    if (@mem && $mem[-1]->{'type'} eq "D" && $l == 1) {
		# increment previous
		$mem[-1]->{'length'} += 1;
	    } else {
		if ($l > 1) {
		    push @mem, {'type' => 'M', 'length' => $l - 1};
		    $covered += $l - 1; 
		}
		push @mem, {'type' => 'D', 'length' => 1};
	    }
	}
    }
}

for (sort keys %reflen) {
    $head .= "\@SQ\tSN:$_\tLN:$reflen{$_}\n";
}

$head .= $pg;

print $head, $body;

# unmapped target
exit if $skip_files;
for my $t (sort keys %$fasta) {
    next if $seen{$t};
    # Query2 4 * 0 0 * * 0 0 TCCTCGCTCTCGACATGTCTCGCACGGCT * AS:i:0 XS:i:0
    print join("\t", $t, 4, "*", 0, 0, "*", "*", 0, 0, "*", $qual) . "\n";    
}

# subroutine
sub read_fasta {
    my ($file) = @_;
    my %hash;
    # Use STDIN if file is '-'
    $file = undef if $file && $file eq '-';
    my $in;
    if ($file && -e $file) {
	open $in, '<', $file || die $!;
    } else {
	$in = *STDIN;
    }
    # Store the sequence id
    my $seqid;
    for (<$in>) {
        # Remove line endings
        s/\R//g;
	# Skip empty lines
	next if /^\s*$/;
	# Check wheter it is an id line
	if (/^>(\S+)/) {
	    $seqid = $1;
	} else {
	    s/\s+//g;
	    # Add to the sequence
	    $hash{$seqid} .= $_;
	}
    }
    close $in;
    return \%hash;
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
