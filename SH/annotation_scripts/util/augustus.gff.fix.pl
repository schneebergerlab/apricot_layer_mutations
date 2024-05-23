#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: augustus.gff.fix.pl
#
#        USAGE: ./augustus.gff.fix.pl  
#
#  DESCRIPTION: a
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 03/27/20 11:13:58
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($in,$out)=@ARGV;
open IN,$in;
open OUT,">$out";
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  my $id = (split /\;/,$t[8])[0];
  $id = (split /=/,$id)[1];
  #if ($t[2] eq "gene") {
    s/ID=/ID=\Q$t[0]\E./g;
    #}else 
    s/Parent=/Parent=\Q$t[0]\E./g;
  print OUT "$_\n";
  
}
close OUT;
close IN;
