#!/usr/bin/env perl
use warnings;
use strict;

use FindBin;                     # locate this script
use lib "$FindBin::RealBin/lib";  # use the lib directory
use biointsam;
use biointbasics;

my $programname = "sam-flip.pl";
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

# Print help if needed
biointbasics::print_help(\@ARGV, $info);

my %ref;
my %query;

open(my $temp, '>', 'temp~') || die $!;

my $header = "true";
my $previousprogram = "";

while(<>) {
    my %hit;
    biointsam::parse_sam($_, \%ref, \%hit); 
    unless (%hit) {
	print {$temp} "$_\n" unless /^\@SQ/;
	if (/^\@PG\tID:(\S+)/) {
	    $previousprogram = $1;
	}
	next;
    }
    if ($header) { # First line after the header section
	my $text = "\@PG\tID:$programname\tPN:$programname";
	$text .= "\tPP:$previousprogram" if $previousprogram;
	print {$temp} $text . "\tVN:$version\tCL:$cmd\n";
	$header = undef;
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
