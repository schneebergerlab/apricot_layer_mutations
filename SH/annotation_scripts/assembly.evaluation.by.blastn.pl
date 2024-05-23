#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: assembly.evaluation.by.2019HR.pl
#
#        USAGE: ./assembly.evaluation.by.2019HR.pl  
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
#      CREATED: 03/02/2020 01:35:36 PM
#     REVISION: ---
#===============================================================================
use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO::blast;
use File::Basename;

my ($assFas,$geneFas,$outdir,$outp)=@ARGV;

##run blastn
system("blastn -query $geneFas -db $assFas -outfmt 6 -out $outdir/$outp.out.m6 -num_threads 20");
system("blastn -query $geneFas -db $assFas -out $outdir/$outp.out -num_threads 20");

open OUT,">$outdir/$outp.summary.out";
my $blastOut = "$outdir/$outp.out";
parseBlast($blastOut,$outdir,$outp);
   
   open IN,"$outdir/$outp.besthit.out";
   my ($total,$total2) = (0,0);   
   while (<IN>) {
   	   chomp;
   	   my @t=split /\t/;
   	   if ($t[1] eq "None") {
   	   	   $total++;	   	  
   	   }elsif ( ($t[8]<95 ) and ($t[9]<85))  {
   	   	   $total2++;	   	  
   	   }   	   
   }
   close IN;
   print OUT "nohit\t$total\n";
   print OUT "partial\t$total2\n";
   close OUT;

sub parseBlast {
  my ($blast_result,$outdir,$o)=@_;
  my $in = new Bio::SearchIO(-format => 'blast', -file => "$blast_result",-best_hit_only=>"TRUE");
  open OUT1,">$outdir/$o.besthit.out";
  open OUT2,">$outdir/$o.besthit.mis.gap.out";
  while( my $result = $in->next_result ) {
    ## $result is a Bio::Search::Result::ResultI compliant object
    my $query=$result->query_name;
    my %hsp;my $score=0;
    my $qleng=$result->query_length();
    my ($subject,$gaps,$mis,$qs,$qe,$ss,$se,$cov,$iden)=("None",0,0,0,0,0,0,0,0);
    my @mis=();my @gaps=();
    my @ins=();my @del=(); ### the subject comapring with final gene seq.
    while( my $hit = $result->next_hit ) {
      while( my $hsp = $hit->next_hsp ) {
      ## $hsp is a Bio::Search::HSP::HSPI compliant object                
        if ( $score < $hsp->score() ) {
          $score=$hsp->score();
          $subject=$hit->name;
          $iden=$hsp->percent_identity;
          $iden=sprintf("%.1f",$iden);
          $qs=$hsp->start("query");
          $qe=$hsp->end("query");
          $ss=$hsp->start("hit");
          $se=$hsp->end("hit");
          $gaps=$hsp->gaps('total');
          $cov=sprintf("%.1f",100*(abs($qe-$qs)+1)/$qleng) ;
          
          @gaps=$hsp->seq_inds('query','gap');
          @ins=@gaps;
          my %gaps;
          foreach (@gaps) {
          	$gaps{$_}=1;
          }
          
          my %mis;          
          @mis=$hsp->seq_inds('query','mismatch');
          foreach (@mis) {
          	$mis{$_}=1;
          }
          
          my @nomatch=$hsp->seq_inds('query','nomatch');
          foreach (@nomatch) {
          	if (!$mis{$_}) {
          	  push @del,$_ if (!$gaps{$_});
          	}
          }
          if (@mis) {
          	$mis=$#mis+1;
          }
          
        }else {
          next;
        }        
      }
    }  
    
    print OUT1 "$query\t$subject\t$mis\t$gaps\t$qs\t$qe\t$ss\t$se\t$cov\t$iden\t$qleng\n";
    
    if (@mis) {
      foreach my $pos (@mis) {
        print OUT2 "$query\t$pos\t$pos\t$subject\t$mis\t$gaps\t$qs\t$qe\t$ss\t$se\t$cov\t$iden\tmis\t$qleng\n";	
      }      
    }
    if (@ins) {
      foreach my $pos (@ins) {
        print OUT2 "$query\t$pos\t$pos\t$subject\t$mis\t$gaps\t$qs\t$qe\t$ss\t$se\t$cov\t$iden\tins\t$qleng\n";	
      }	
    }
    if (@del) {
      foreach my $pos (@del) {
        print OUT2 "$query\t$pos\t$pos\t$subject\t$mis\t$gaps\t$qs\t$qe\t$ss\t$se\t$cov\t$iden\tdel\t$qleng\n";	
      }	
    }    
  }
  close OUT1;close OUT2;
}
