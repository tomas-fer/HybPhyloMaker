#!/lusr/bin/perl -w
#https://github.com/smirarab/global/blob/master/src/perl/convert_unit_to_100_support.pl

use strict;
use warnings;

#input: tree with support values in range 0-1
#output: tree with support values in range 0-100

my $tree = $ARGV[0];

sub printUsage {
    my $msg = "$0 <tree>\n";
    die $msg;
}

printUsage() unless(defined($tree));

my $output = "$tree.100";
my $tree_contents = `cat $tree`;
$tree_contents =~ s/\)(\d\.?\d*)/scale_br($1);/ge;
open(OUT, ">", $output) or die "can't open $output: $!";
 print OUT $tree_contents;
close(OUT);

print "output at $output\n";


sub scale_br {
 my ($br) = @_;
 return ")".(int($br*100));
}