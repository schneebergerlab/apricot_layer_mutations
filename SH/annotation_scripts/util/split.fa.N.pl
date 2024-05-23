#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: split.fa.N.pl
#
#        USAGE: ./split.fa.N.pl  
#
#  DESCRIPTION: split a fasta file into multiple smaller files
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 03/18/20 17:53:56
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;

my ($in,$outdir,$p,$num)=@ARGV;

my $seqin = new Bio::SeqIO(-format=>"fasta",-file=>$in);
my $n=0;my $k = 1;
open OUT,">$outdir/$p.$k.fa";
while (my $seqobj = $seqin->next_seq() ){
   my $id = $seqobj->id();
   my $seq = $seqobj->seq();
   if ($n==$num) {
     close OUT;
     $k+=1;
     $n=1;
     open OUT,">$outdir/$p.$k.fa";
     print OUT ">$id\n$seq\n";
   }else {
     $n+=1;
     print OUT ">$id\n$seq\n";
   }
}
close OUT;
