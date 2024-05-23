#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: orthofinder.blastp.id.replace.pl
#
#        USAGE: ./orthofinder.blastp.id.replace.pl  
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
#      CREATED: 03/17/20 11:16:59
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($idFile,$blastp,$out)=@ARGV;


open IN,$idFile;
my %id;
while (<IN>) {
  chomp;
  s/\://g;
  my @t = (split / /);
  #print $t[0];exit;
  $id{$t[0]}=$t[1];
}
close IN;

open IN,$blastp;
open OUT,">$out";
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  $t[0] = $id{$t[0]};
  $t[1] = $id{$t[1]};
  my $line = join("\t",@t);
  print OUT "$line\n";
}
close IN; close OUT;



