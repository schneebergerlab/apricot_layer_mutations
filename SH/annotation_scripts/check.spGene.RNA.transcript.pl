#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: check.spGene.RNA.transcript.pl
#
#        USAGE: ./check.spGene.RNA.transcript.pl  
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
#      CREATED: 04/04/20 11:18:58
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


my ($spBed, $inGff, $rnaBed, $outdir,$o) = @ARGV;


open IN,$spBed;
my %spGenes;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  $spGenes{$t[4]}=1;
}
close IN;

open IN,$inGff;
open OUT,">$outdir/$o.spGene.exon.bed";
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  if ($t[2] eq "exon") {
    my $pa = (split /\;/,$t[8])[-1];
     $pa = (split /\=/,$pa)[1];
    my $id = (split /\./,$pa)[0];
    if ($spGenes{$id}) {
      print OUT "$t[0]\t$t[3]\t$t[4]\t$t[6]\t$id\n";
      $spGenes{$id}+=$t[4]-$t[3];
    }
  }
}
close IN;close OUT;


open IN,"intersectBed -a $outdir/$o.spGene.exon.bed -b $rnaBed -wo |";
my %olap;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  $olap{$t[4]}+=$t[-1];
}
close IN;

open OUT,">$outdir/$o.spGene.intersect.RNAseq.transcript.stats";
foreach my $g (sort keys %spGenes) {
  if ($olap{$g}) {
    my $p = $olap{$g}/$spGenes{$g};
    print OUT "$g\t$p\n";
  }else {
     print OUT "$g\t0\n";
  }
}

close OUT;

