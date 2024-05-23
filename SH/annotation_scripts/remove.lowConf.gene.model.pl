#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: remove.lowConf.gene.model.pl
#
#        USAGE: ./remove.lowConf.gene.model.pl  
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
#      CREATED: 03/26/20 10:51:51
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($canAdded,$ungrBed, $athOrt,$rnaGff,$annGff,$outGff,$rmOut) = @ARGV;

# geth gene ungrouped genes which could be futher added (add.false.spGenes.pl)
open IN,$canAdded;
my %add;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  $add{$t[0]} = 1;
}
close IN;


open IN,$ungrBed;
my %ungr;
my $p = "";
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  $ungr{$t[4]} = 1;
  $p = substr($t[4],0,3);
}
close IN;

open IN,$athOrt;
my %ort;
while (<IN>) {
  chomp;
  my @t = (split / /);
  shift @t;
  my @genes ; my @ath;
  foreach my $g (@t) {
    if (substr($g,0,3) eq $p) {
       push @genes, $g;
    }elsif (substr($g,0,2) eq "AT") {
       push @ath, $g;
    }
  }
  next if $#ath < 0; 
  foreach my $g (@genes) {
  	$g = (split /\./, $g)[0];
    $ort{$g}=1;
  }
}
close IN;

open IN,"intersectBed -a $ungrBed -b $rnaGff -wao |";
print "intersectBed -a $ungrBed -b $rnaGff -wao\n";
my %rna;
while (<IN>) {
   chomp;
   my @t = (split /\t/);
   next if ($t[-1]==0);
   $rna{$t[4]}+= $t[-1];
}
close IN;

open IN,$annGff;
open OUT,">$outGff";
open OUT1,">$rmOut";
my $flag = 0;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  if ($t[2] eq "gene") {
    my $id = (split /\;/,$t[8])[0];
    $id = (split /\=/,$id)[1];
    if ( ! defined $ungr{$id} or $add{$id}) {
       $flag = 1;
       print OUT "$_\n";
    }else {
      if (!$ort{$id} and !$rna{$id}) {
        $flag = 0;
        print OUT1 "$id\n";
      }elsif (!$ort{$id} and $rna{$id}< 250) {
        $flag = 0;
        print OUT1 "$id\n";
      }else {
         $flag = 1;
         print OUT "$_\n";
      } 
    }
  }else {
    if ($flag == 1) {
      print OUT "$_\n";
    }
  }
}
close IN;close OUT;close OUT1;


