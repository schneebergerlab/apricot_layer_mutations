#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: split.genome.for.infernal.pl
#
#        USAGE: ./split.genome.for.infernal.pl  
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
#      CREATED: 08/30/2017 01:02:51 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;

my ($in,$out)=@ARGV;

my $seqin=new Bio::SeqIO(-format=>"fasta",-file=>"$in");
my $seqout=new Bio::SeqIO(-format=>"fasta",-file=>">$out");
while (my $seqobj=$seqin->next_seq()) {
  my $id=$seqobj->id();
  my $len=$seqobj->length();
  my $seq=$seqobj->seq();
  if ($len>10000000) {
    my $n=int($len/10000000);
    foreach my $i (0..$n-1) {
      my $s=substr($seq,$i*10000000,10000000);
      my $seqo=new Bio::Seq(-id=>"$id\_$i",-seq=>"$s");
      $seqout->write_seq($seqo);
    }
    my $s=substr($seq,$n*10000000,$len-$n*10000000);
    my $seqo=new Bio::Seq(-id=>"$id\_$n",-seq=>"$s");
    $seqout->write_seq($seqo);
    
  }else {
    $seqout->write_seq($seqobj);
  }
}


