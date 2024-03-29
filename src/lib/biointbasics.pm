#!/usr/bin/perl

use strict;
use FindBin;                     # locate this script
use lib "$FindBin::RealBin/../lib";  # use the parent directory


# Biont basics: frequently used basic perl solutions
package biointbasics;

sub print_help {
    # Print out the usage to STDERR
    # Takes in the @ARGV as input
    my ($args, $info, $error) = @_;
    $error = "" unless $error;
    my $help_message = join("\n", map{ $info->{$_} } qw/usage description options/);
    if ($error) {
	# Die on errors and add help message
	die $error . $help_message;
    }
    for (@$args) {
	if (/^-?-h(elp)?$/) {
	    # Print help to STDOUT and exit properly
	    print $help_message;
	    exit;
	}
    }
}

sub contains {
    # Check if all the elements of an array are present in a hash
    my ($list, $ref) = @_;
    return scalar(@$list) == scalar( &overlap($list, $ref) );
}

sub overlap {
    # Check which elements of an array are present in a hash
    my ($list, $ref) = @_;
    my @found;
    for (@$list) {
	push @found, $_ if $ref->{$_};
    }
    return @found;
}

sub read_fasta {
    # Convert FASTA string into a hash with IDs for keys and sequences
    #  as values and stores the original order in an array
    # This subroutine requires three arguments:
	#	1) filehandle for the FASTA file
	#	2) a hash reference to store the sequences in
	#	3) an array reference to store the IDs in the same
	#          order as the original file 
    # If an ID line is present multiple times then a warning is printed
    #  to STDERR
    my ($hash, $list, $file) = @_;
    # Use STDIN if file is '-'
    $file = undef if $file && $file eq '-';
    my $in;
    if ($file && -e $file) {
	open $in, '<', $file || die $!;
    } else {
	$in = *STDIN;
    }
    # Store the sequence id
    my $seq_id;
    for (<$in>) {
        # Remove line endings
        s/\R//g;
	# Skip empty lines
	next if /^\s*$/;
	# Check wheter it is an id line
	if (/>(.*)/) {
	    # Save the id and the definition and store it in the array
	    $seq_id = $1;
	    print {*STDERR} "WARNING: <$seq_id> is present in multiple copies\n" if $hash->{$seq_id};
	    push @$list, $seq_id;
	} else {
	    # If there was no id lines before this then throw an error
	    unless (defined $seq_id) {
		print "Format error in FASTA file! Check the file!\n";
		last;
	    }
	    # Remove white space
	    s/\s+//g;
	    # Add to the sequence
	    $hash->{$seq_id} .= $_;
	}
    }
    close $in;
}

1;
