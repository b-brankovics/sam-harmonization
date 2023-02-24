#!/usr/bin/perl

use strict;
use FindBin;                     # locate this script
use lib "$FindBin::RealBin/../lib";  # use the parent directory


# Biont SAM package
package biointsam;

sub score {
    # update alignment score
    my ($hit) = @_;
    # Get alginment data
    my $aln = biointsam::parse_cigar($hit->{'CIGAR'}, $hit->{'FLAG'}, $hit->{'SEQ'});

    my $score = 0;
    # Scoring rules
    my $si = 5;  # identical
    my $sm = -4; # mismatch or gapextend
    my $sg = -8; # gapopen (-12) - gapextend (-4)
    my $gapcount = () = $hit->{'CIGAR'} =~ /[ID]/g;
    # exonerate score = 5 * iden - 4 * NM:i - 8 * gapcount 
    #                 = 5 * ($aln->{'length'} - NM:i) - 4 * NM:i - 8 * gapcount 
    #                 = 5 * $aln->{'length'} - 9 * NM:i - 8 * gapcount
    my $exoscore = 5 * $aln->{'length'} - 9 * $hit->{'NM:i'} - 8 * $gapcount;
    my $identical = $aln->{'length'} - $hit->{'NM:i'};
    my $in = $aln->{'insertion'};
    my $del = $aln->{'deletion'};
    my $gappos = $in + $del;
    my $mismatch = $hit->{'NM:i'} - $gappos;
    #$score = $identical - 2 * $mismatch - 2 * ($in + $del) - $gapcount;
    unless ($hit->{'AS:i'} && $hit->{'AS:i'} == $exoscore) {
	$hit->{'AS:i'} = 'NA' unless $hit->{'AS:i'};
	print STDERR "Updating Alignment score from " . $hit->{'AS:i'} . " to " . $exoscore ."\n";
    }
    $hit->{'AS:i'} = $exoscore;
    return;
}

sub flip {
    # Takes a ref to a hit hash and exchanges Ref with Query
    my ($hash, $ref) = @_;
    # Get alignment info
    my $aln = biointsam::parse_cigar($hash->{'CIGAR'}, $hash->{'FLAG'}, $hash->{'SEQ'});
    my $unmapped1 = $hash->{'POS'} - 1;
    my $mapped = $aln->{'length'} - $aln->{'insertion'};
    my $unmapped2 = $ref->{$hash->{'RNAME'}} - $mapped - $unmapped1;
    my @clips = ($unmapped1, $unmapped2);
    
    $hash->{'POS'} = $aln->{'start'};
    # Flip names
    my $name = $hash->{'QNAME'};
    $hash->{'QNAME'} = $hash->{'RNAME'};
    $hash->{'RNAME'} = $name;

    my $old = $hash->{'CIGAR'};
    die "ERROR: unexpected CIGAR element '$1' in $old for $name\n" if $old =~ /(\d*[^0-9SHDIM=X])/;
    $hash->{'CIGAR'} = "";
#    print STDERR "$mapped @clips $old ";
    $old =~ s/^\d+[SH]//;
    $old =~ s/\d+[SH]$//;
    if ($hash->{'FLAG'} & 16) {
	@clips = reverse @clips;
	$old =~ s/(\d+\D)/$1 /g;
	my $new = join("", reverse split/ /, $old);
	$old = $new;
    }
    if ($clips[0]) {
	$hash->{'CIGAR'} .= $clips[0] . "S";
    }
    $old =~ tr/ID/DI/;
    $hash->{'CIGAR'} .= $old;
    if ($clips[1]) {
	$hash->{'CIGAR'} .= $clips[1] . "S";
    }
#    print STDERR $hash->{'CIGAR'} . "\n";
    # Cannot flip the sequence unless we can read the ref file
    $hash->{'SEQ'} = "*";
}
    
sub sam_string {

    my ($hash) = @_;
    my @field = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;    
    my %done;
    my @list;
    for (@field) {
	push @list, $hash->{$_};
	$done{$_}++;
    }
    # Bring NM:i and AS:i to the front if defined
    for ('NM:i', 'AS:i') {
	next unless defined $hash->{$_};
	push @list, $_ . ":" . $hash->{$_};
	$done{$_}++;
    }
    for (sort keys %$hash) {
	next if $done{$_};
	push @list, $_ . ":" . $hash->{$_};
    }
    return join("\t", @list);
}

sub parse_sam {
    # Reads a SAM file line by line and converts each entry to a hash
    my ($line, $ref, $hash) = @_;
    $_ = $line;
    my @field = qw/QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL/;    
    s/\R+//;
    my @col = split/\t/;
    if (/^\@/) { # header
	my %h;
	for (@col[1..$#col]) {
	    if (/^([^:]+):(.*)/) {
		$h{$1} = $2;
	    } else {
		print STDERR "Warning: could not parse header element '$_'\n";
	    }
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

sub parse_cigar {
    my ($cigar, $bitflag, $seq) = @_;
    $seq = undef if $seq eq '*';
    my $len = 0;
    my $insert = 0;
    my $deletion = 0;
    my $start = 1;
    my $unmapped = 0;
    # How long are the soft clipped regions
    my $soft5 = 0;
    my $soft3 = 0;
    # Remove clipped ends and store info
    if ($cigar =~ s/^(\d+)([SH])//) {
	my $length = $1;
	$start += $length unless $bitflag & 16;
	$unmapped += $length;
	if ($2 eq 'S' && $seq) {
	    #$seq =~ s/^.{$length}//;
	    $seq = substr($seq, $length);
	    $soft5 = $length;
	}
    }
    if ($cigar =~ s/(\d+)([SH])$//) {
	my $length = $1;
	$start += $length if $bitflag & 16;
	$unmapped += $length;
	if ($2 eq 'S' && $seq) {
	    #$seq =~ s/.{$length}$//;
	    $seq = substr($seq, 0, -$length);
	    $soft3 = $length;
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
    return {"start" => $start, # first position of query mapped
	    "length" => $len, # alignment length
	    "insertion" => $insert,
	    "deletion" => $deletion,
	    "unmapped" => $unmapped, # total unmapeed positions
	    "reverse" => $bitflag & 16, # 0 or 16
	    "seq" => $seq, # sequence that is mapped
	    "soft5" => $soft5, # sequence soft clipped at the start
	    "soft3" => $soft3, # sequence soft clipped at the end
    };
}

sub parse_blast {
    # Parse the BLAST output
    my ($output, $cmd, $version) = @_;
    my %reflen;
    my $head = "";
    my $body = "";
    my $mapq = 255;
    my $unmapped = "";

    my $mem;
    for (split/\R+/, $output) {
	if (/^#/) {
	    if (/^#\s+Query:\s+(\S+)/) {
		$mem = $1;
	    } elsif (/^#\s+0\s+hits?\s+found/) {
		# No hits
		$unmapped .= $mem . "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
	    }
	    next;
	}
	#    print "$_\n";
	#    next;
	#    s/\R+//;
	# sseqid slen sstart send qseqid qlen qstart qend sstrand mismatch gaps btop sseq qseq
	my ($r, $rlen, $rab, $rae, $t, $tlen, $tab, $tae, $ts, $mismatch, $gaps, $btop, $rseq, $tseq, $score, $bitscore, $evalue) = split/\t/;
	# Set refID to match the ID in the FASTA file
	if ($r =~ /\|([^|]+)\|$/) {
	    $r = $1;
	}
	$reflen{$r} = $rlen;
	my $bitflag = 0;
	my $seq = $tseq;
	# Remove gaps
	$seq =~ s/-//g;
	if ($ts eq "minus") {
	    # Reverse match
	    $bitflag += 16;
	    $seq = &reverse_seq($seq);
	}
	# Sort positions
	($rab, $rae) = sort{$a<=>$b} ($rab, $rae);
	($tab, $tae) = sort{$a<=>$b} ($tab, $tae);
	my $qual = "*";
	#
	my $edit = $mismatch + $gaps;
	
	my $pos = $rab;
	my $clip1 = $tab - 1;
	my $clip2 = $tlen - $tae;
	my $start = "";
	my $end = "";
	$start = $clip1 . "H" if $clip1;
	$end = $clip2 . "H" if $clip2;
	($start, $end) = ($end, $start) if $bitflag & 16;
	my $cigar = $start;
	#$exocigar =~ s/\s//g;
	my @mem;
	while ($btop) {
	    $btop =~ s/^((\D{2})|(\d+))//;
	    my $x = $1;
	    unless ($x) {
		print STDERR "problem with BTOP: '$_'\n";
		last;
	    }
	    my $type;
	    my $l;
	    if ($x =~ /^\d+$/) {
		$type = "M";
		$l = $x;
	    } else {
		$l = 1;
		# X->query Y->subject/ref
		my $y = substr($x, 1, 1, '');
		if ($x eq "-") {
		    # missing from query => D
		    $type = "D";
		} elsif ($y eq "-") {
		    # missing from subject/ref => I
		    $type = "I";
		} else {
		    # mismatch
		    $type = "M";
		}
	    }
	    if (@mem && $mem[-1]->{'type'} eq $type) {
		$mem[-1]->{'length'} += $l;
	    } else {
		push @mem, {'type' => $type, 'length' => $l};
	    }
	    
	}
	while (@mem) {
	    $_ = shift @mem;
	    $cigar .= $_->{'length'} . $_->{'type'};
	}
	# Add remaining match
	
	$cigar .= $end;
	
	$body .= join("\t", $t, $bitflag, $r, $pos, $mapq, $cigar, "*", 0, 0, $seq, $qual, "NM:i:" . $edit, "AS:i:$score", "BS:f:$bitscore", "EV:f:$evalue") . "\n";
	
    }
    for (sort keys %reflen) {
	$head .= "\@SQ\tSN:$_\tLN:$reflen{$_}\n";
    }
    
    $head .= "\@PG\tID:BLAST\tPN:blastn\tVN:$version\tCL:$cmd\n";
    
    print $head, $body, $unmapped;

}

sub parse_aln {
    # Takes two sequences that are aligned and returns the CIGAR and edit distance
    my ($seq1, $seq2) = @_;
    my $total = length $seq1;

    die "ERROR: sequences have different lengths\n" unless $total == length($seq2);
    my %hit;
    my $edit = 0;
    my $type = "S";

    $hit{"QNAME"} = "*";
    $hit{"FLAG"} = 0;
    $hit{"RNAME"} = "*";
    $hit{"POS"} = 1;
    $hit{"MAPQ"} = 255;
    # CIGAR comes later
    $hit{"RNEXT"} = '*';
    $hit{"PNEXT"} = 0;
    $hit{"TLEN"} = 0;
    $hit{"SEQ"} = $seq2;
    $hit{"SEQ"} =~ s/-//g;
    $hit{"QUAL"} = '*';

    my $len;
    my $cigar = "";
    # Loop through each position in the alignment
    for my $p (1..$total) {
	my $base1 = substr($seq1, $p - 1, 1);
	my $base2 = substr($seq2, $p - 1, 1);
	if ($type eq "S") {
	    if ($base1 eq "-" && $hit{"POS"} == 1) {
		# Ref hasn't started yet => Soft clip (S)
		&cigar_update(\$type, \$len, "S", \$cigar);
		next; # skip adding an edit distance
	    } elsif ($base1 ne "-" && $base2 eq "-") {
		# Ref started, but Q has gap => starting position
		$hit{"POS"}++;
		next; # skip adding an edit distance
	    } else {
		$cigar .= $len . $type if $len;
		$len = 0;
		$type = "";
	    }
	}
	if ($base1 eq $base2) {
	    next if $base1 eq '-';
	    # Match (M)
	    &cigar_update(\$type, \$len, "M", \$cigar);
	    next; # skip adding an edit distance
	} elsif ($base1 eq "-") {
	    # Insertion (I)
	    &cigar_update(\$type, \$len, "I", \$cigar);
	} elsif ($base2 eq "-") {
	    # Deletion (D)
	    &cigar_update(\$type, \$len, "D", \$cigar);
	} else {
	    # Mismatch (M)
	    &cigar_update(\$type, \$len, "M", \$cigar);
	}
	$edit++;
    }
    if ($type eq "I" && substr($seq1, -1, 1) eq "-") {
	# Ref already ended => Soft clip (S)
	$type = "S";
	$edit -= $len;
    } elsif ($type eq "D") {
	$type = "";
	$edit -= $len;
    }
    # Update the final stretch
    &cigar_update(\$type, \$len, "", \$cigar);
    $hit{"CIGAR"} = $cigar;
    $hit{"NM:i"} = $edit;

    return \%hit;
}

sub cigar_update {
    # Update cigar whether type is switched
    # Increment alignment by one
    #  - type: the current CIGAR tag M, I, D, etc.
    #  - len: length of the current elemnet (e.g. 5 in 5M)
    #  - new: the CIGAR tag of the latest nucleotide
    #  - cigar: the CIGAR string so far
    my ($type, $len, $new, $cigar) = @_;
    if ($$type && $$type eq $new) {
	$$len++;
    } else {
	$$cigar .= $$len . $$type if $$type;
	$$type = $new;
	$$len = 1;
    }
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


# sub new {
#     my ($class, %args) = @_;
#     return bless \%args, $class;
# }

# sub get_provenance {
# 
# }

1;
