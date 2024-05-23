#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: dat2fasta.pl
#
#        USAGE: ./dat2fasta.pl  
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
#      CREATED: 01/28/2020 01:41:25 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($in,$out)=@ARGV;

open IN,$in;
open OUT,">$out";
my $id=""; my $seq="";
while (<IN>) {
  chomp;
  if (/^ID/) {
    s/\s+/\t/g;
    my @t = (split /\t/,$_);
    $id = $t[1];
  }elsif (/^SQ/) {
    while (<IN>) {
      chomp;
      if (/^\//) {
        print OUT ">$id\n$seq\n";
        $id="";$seq="";
        last;
      }
      s/\s+//g;
      $seq.=$_;
    }
  }
}
close IN;close OUT;
