#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: gene.model.validate.pl
#
#        USAGE: ./gene.model.validate.pl  
#
#  DESCRIPTION: 
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 03/25/20 15:07:47
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;

my ($gff,$geneFa,$protFa,$cdsFa,$out) = @ARGV;

open IN,$gff;
my $id = ""; my %len;
my %gene1; my %gff;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  if ($t[8]=~/trans/) {
  	  next;	
  	}
  if ($t[2] eq "gene") {
     $id = (split /\;/,$t[8])[0];
    $id = (split /\=/,$id)[1];
    #print $id, "\n";exit;
    $gene1{$id}=$t[4]-$t[3]+1;
  }
  elsif ($t[2] eq "CDS"){
    $len{$id}+=$t[4]-$t[3] + 1;
  }
}
close IN;

my $seqin = new Bio::SeqIO(-format=>"fasta",-file=>"$geneFa");
my %gene2; my %dif ;
while (my $seqobj = $seqin->next_seq()) {
  my $id = $seqobj->id();
  $id = (split /\./,$id)[0];
  my $len = $seqobj->length();
  $gene2{$id}=$len;
  if ($len != $gene1{$id}) {
     print "diff length\t$id\t$gene1{$id}\t$len\n";
     $dif{$id}=1;
  }
}

$seqin = new Bio::SeqIO(-format=>"fasta",-file=>"$cdsFa");
my %st1;
while (my $seqobj=$seqin->next_seq()) {
  my $id = $seqobj->id();
  $id = (split /\./,$id)[0];
  my $seq = $seqobj->seq();
  if ( substr($seq,length($seq)-3,3) eq "TGA" or 
  substr($seq,length($seq)-3,3) eq "TAA" or 
  substr($seq,length($seq)-3,3) eq "TAG"  
  )  {
    next;
  }else {
    $st1{$id}=1;
    print "not stop-codon\t$id\n";
  }
}

$seqin = new Bio::SeqIO(-format=>"fasta",-file=>"$protFa");
my %st2;
while (my $seqobj = $seqin->next_seq()) {
  my $id = $seqobj->id();
  $id = (split /\./,$id)[0];
  my $seq = $seqobj->seq();
  #if ($seq!~/\.$/) {
  #  print "no stop\t$id\n";
  #  $st1{$id}=1;
  #}
  my $s = substr($seq,0,length($seq)-1);
  if ($s=~/\*/) {
    print "inside stop\t$id\n";
    $st2{$id}=1;
  }
}

open IN,$gff;
open OUT,">$out";
my $flag;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  if ($t[2] eq "gene") {
  	if ($t[8]=~/trans/) {
  	  $flag = 1;
  	   print OUT "$_\n";
  	  next;	
  	}
    my $id = (split /\;/,$t[8])[0];
    $id = (split /\=/,$id)[1];
    if ($len{$id} % 3 != 0) {
       print "no 3 times length of CDS\t$id\n";
       $flag = 0;
    }elsif ($dif{$id} or $st1{$id} or $st2{$id}) {
       $flag = 0;
    }else {
       $flag = 1;
       print OUT "$_\n";
    }
  }else {
     if ($flag == 1) {
       print OUT "$_\n";
     }
  }
}
close IN;close OUT;

