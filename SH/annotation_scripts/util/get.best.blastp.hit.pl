#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: get.best.blastp.hit.pl
#
#        USAGE: ./get.best.blastp.hit.pl  
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
#      CREATED: 09/13/2018 05:43:24 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;

my ($blastp,$seq1,$seq2,$gff,$ann,$out)=@ARGV;

my $leng1=getLeng($seq1);
my $leng2=getLeng($seq2);
my %leng1=%{$leng1};
my %leng2=%{$leng2};
open IN,$ann;
my %ann;
while (<IN>) {
  chomp;
  my @t=(split /\t/);
  $ann{$t[0]}=$t[1];
}
close IN;

open IN,$blastp;
my %blast;
my %iden;
while (<IN>) {
  chomp;
  my @t=(split /\t/);
  $t[0]=(split /\./,$t[0])[0];
  $t[1]=(split /\./,$t[1])[0];
  my $cov1=($t[7]-$t[6])/$leng1{$t[0]};
  my $cov2=abs($t[9]-$t[8])/$leng2{$t[1]};
  if ($blast{$t[0]}{$t[1]}) {
    $blast{$t[0]}{$t[1]}+=$t[2]*$cov1*$cov2;
    if ($iden{$t[0]}{$t[1]} < $t[2]) {
      $iden{$t[0]}{$t[1]}=$t[2];
    }
  }else {
    $blast{$t[0]}{$t[1]}=$t[2]*$cov1*$cov2;
    $iden{$t[0]}{$t[1]}=$t[2];
  }
  
}
close IN;

open IN,$gff;
open OUT,">$out";
while (<IN>) {
  chomp;
  my @t=(split /\t/);
  if ($t[2] eq "gene") {
    my $id=(split /\;/,$t[8])[0];
    $id=(split /\=/,$id)[1];
    if ($blast{$id}) {
      my %tt = %{$blast{$id}};
      my @tmp = sort {$tt{$b}<=>$tt{$a}} keys %tt;
      my $hit = $tmp[0];
      my $iden =$iden{$id}{$hit};
      $t[8]=$t[8].";AtBestHit=$hit;Identity=$iden;Ann=$ann{$hit}";
    }else {
      $t[8]=$t[8].";AtBestHit=None";
    }
  }
  my $line=join("\t",@t);
  print OUT "$line\n";
}
close IN;close OUT;


sub getLeng {
  my ($file)=@_;
  my $seqin=new Bio::SeqIO(-format=>"fasta",-file=>"$file");
  my %leng;
  while (my $seqobj=$seqin->next_seq()) {
    my $id=$seqobj->id();
    $id=(split /\./,$id)[0];
    my $leng=$seqobj->length();
    $leng{$id}=$leng;
  }
  return \%leng;
}

