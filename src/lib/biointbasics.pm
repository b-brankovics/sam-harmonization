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


1;
