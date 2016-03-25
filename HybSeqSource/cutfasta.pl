#!/usr/bin/env perl

# ******************************************************************************************************
# * Tomas Fer, Dept. of Botany, Charles University, Prague, Czech Republic, 2015                       *
# * tomas.fer@natur.cuni.cz                                                                            *
# * based on: http://stackoverflow.com/questions/16553004/extracting-fasta-sequences-based-on-position *
# ******************************************************************************************************

# Cuts specified range from (aligned) multiple fasta file (e.g., bases from 100 to 200)
# Usage perl cutfasta.pl input.fas from to > output.fas, e.g. perl cutfasta.pl input.fas 100 200
use warnings;
use strict;

my ($adn, $l, $header);
my $end   = pop;
my $start = pop;
while ( <> ) { 
    chomp;

    ## First line is known, a header, so print it and process next one.
    if ( $. == 1 ) { 
		printf qq|%s_%s\n|, $_, $start . "_" . $end;
        next;
    }   

    ## Concat adn while not found a header.
    if ( '>' ne substr $_, 0, 1 ) { 
        if ( ! $l ) { $l = length }
        $adn .= $_; 
        if ( ! eof ) { next }
    }   
    else {
        $header = sprintf qq|%s_%s\n|, $_, $start . "_" . $end;
    }   

    ## Extract range from-to and insert newlines to set same length of 
    ## line than before.
    my $s = substr $adn, $start, $end - $start;
    $s =~ s/(.{$l})/$1\n/g;
    printf qq|%s\n|, $s; 
    undef $adn;

    ## If not end of file, print the header of the following adn.
    if ( ! eof ) { print $header }
}