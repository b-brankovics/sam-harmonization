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
    "A wrapper for exonerate that produces valid SAM files as output.\n" .
    "\tThe affine:local model (gapped alignment) will be used.\n";
my $usage = 
    "Usage:\n\t$0 [OPTIONS] <reference> <query>\n";
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

my ($ref, $query) = @ARGV;

for my $file ($ref, $query) {
    if (! defined $file) {
	biointbasics::print_help(\@ARGV, $info,
				 "ERROR: You need to specify exactly two input FASTA files.\n");
    } elsif (! -e $file || -z $file) {
	biointbasics::print_help(\@ARGV, $info,
				 "ERROR: '$file' does not exist or it is empty.\n");
    }
}

my $versiontext = `exonerate -v`;
$versiontext =~ /version (\S+)/;
my $version = $1;
my $cmd = "exonerate $ref $query -m affine:local --verbose 0 --showalignment no --showvulgar no --ryo ";
#$cmd .= '"%qi\t%ql\t%qab\t%qae\t%qS\t%ti\t%tl\t%tab\t%tae\t%tS\t%C\t%em\n"';
$cmd .= '"%qi\t%ql\t%qab\t%qae\t%qS\t%ti\t%tl\t%tab\t%tae\t%tS\t%C\t%em\t%s\t%tas\n"';

my %reflen;
my $head = "\@PG\tID:exonerate\tPN:exonerate\tVN:$version\tCL:$cmd\n";
my $body = "";

my $mapq = 255;

my $output = `$cmd`;
$output =~ s/\n([A-Za-z]+)/$1/g;
#open my $input, "$cmd |";
for (split/\R+/, $output) {
#    print "$_\n";
#    next;
#    s/\R+//;
    my ($r, $rlen, $rab, $rae, $rs, $t, $tlen, $tab, $tae, $ts, $exocigar, $mismatch, $score, $seq) = split/\t/;
    $reflen{$r} = $rlen;
    my $bitflag = 0;
    if ($rs ne $ts) {
	# Reverse match
	$bitflag += 16;
	if ($ts eq "-") {
	    $rab++;
	    $tae++;
	} else {
	    $rae++;
	    $tab++;
	}
    } elsif ($rs eq "-") {
	# Both are reversed
	$rae++;
	$tae++;
    } else {
	# Not reversed
	$rab++;
	$tab++;
    }
    # Sort positions
    ($rab, $rae) = sort{$a<=>$b} ($rab, $rae);
    ($tab, $tae) = sort{$a<=>$b} ($tab, $tae);
    my $qual = "*";
    $seq = "" unless $seq;
    my $edit = $mismatch;
	
    if ( ($rab > $rae && $tab < $tae) || ($rab < $rae && $tab > $tae) ) {
	
    }
    my $pos = $rab;
    my $clip1 = $tab - 1;
    my $clip2 = $tlen - $tae;
    my $start = "";
    my $end = "";
    $start = $clip1 . "H" if $clip1;
    $end = $clip2 . "H" if $clip2;
    ($start, $end) = ($end, $start) if $bitflag & 16;
    my $cigar = $start;
    $exocigar =~ s/\s//g;
    while ($exocigar) {
	if ($exocigar =~ s/^(\D)(\d+)//) {
	    my ($type, $length) = ($1, $2);
	    # Correct exonerate on I and D
	    $type =~ tr/DI/ID/;
	    $cigar .= $length . $type;
	    if ($type =~ /^[DIN]$/) {
		$edit += $length;
	    }
	} else {
	    print STDERR "problem with '$exocigar'\n";
	    last;
	}
    }
    $cigar .= $end;

    $body .= join("\t", $t, $bitflag, $r, $pos, $mapq, $cigar, "*", 0, 0, $seq, $qual, "NM:i:" . $edit, "AS:i:" . $score) . "\n";
    
}
for (sort keys %reflen) {
    $head .= "\@SQ\tSN:$_\tLN:$reflen{$_}\n";
}

print $head, $body;
