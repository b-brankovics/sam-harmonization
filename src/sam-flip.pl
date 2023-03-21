#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;

my %ref;
my %query;

open(my $temp, '>', 'temp~') || die $!;


while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    unless (%hit) {
	print {$temp} "$_\n" unless /^\@SQ/;
	next;
    }
    next if $hit{"FLAG"} & 4 || $hit{'RNAME'} eq "*";
    # Get alginment data
    my $aln = biointsam::parse_cigar($hit{'CIGAR'}, $hit{'FLAG'}, $hit{'SEQ'});
    $query{ $hit{'QNAME'} } = $aln->{'length'} - $aln->{'deletion'} + $aln->{'unmapped'};
    biointsam::flip(\%hit, \%ref);
    print {$temp} biointsam::sam_string(\%hit), "\n";
}
#print "\n";
close $temp;
for (sort{$a cmp $b} keys %query) {
    print  join("\t", '@SQ', "SN:$_", "LN:$query{$_}"), "\n";
}

#print "\n";
open(my $temp2, '<', 'temp~') || die $!;
while (<$temp2>) {
    print;
}
close $temp2;
unlink 'temp~';
