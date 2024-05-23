#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: srt.gff.pl
#
#        USAGE: ./srt.gff.pl  
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
#      CREATED: 09/25/2017 02:54:25 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;

my ($gff,$out)=@ARGV;

open IN,$gff;
my %genes1;my %genes2;
my @tmp ;
my $start=0;
my $id="";
my $ch="";
while (<IN>) {
  chomp;
  my @t=split /\t/;
  if ($t[2] eq "gene") {
    if ($id ne "") {
      $genes1{$ch}{$id}=$start;
      $genes2{$ch}{$id}=[@tmp];
      $ch = $t[0];
      $id = (split /\;/,$t[8])[0];
      $id = (split /\=/,$id)[1];
      @tmp = ();
      push @tmp,$_;
      $start = $t[3];
    }else {
      $start = $t[3];
      $id = (split /\;/,$t[8])[0];
      $id = (split /\=/,$id)[1];
      push @tmp,$_;
      $ch = $t[0];
    }
  }elsif (($t[2] eq "CDS") or ($t[2] eq "exon") or ($t[2] eq "mRNA")) {
    push @tmp,$_;
  }
}
$genes1{$ch}{$id}=$start;
$genes2{$ch}{$id}=[@tmp];

open OUT,">$out";
foreach my $chr (sort keys %genes1) {
  foreach my $id (sort {$genes1{$chr}{$a}<=>$genes1{$chr}{$b}} keys %{$genes1{$chr}}) {
    if (!$genes2{$chr}{$id}) {
      print  "$chr\t$id\n";
      exit;
    }
    my $line=join("\n",@{$genes2{$chr}{$id}});
    print OUT "$line\n";
  }
}
close OUT;
