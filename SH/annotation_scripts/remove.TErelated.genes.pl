#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: remove.TErelated.genes.pl
#
#        USAGE: ./remove.TErelated.genes.pl  
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
#      CREATED: 08/21/2017 05:20:47 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO::blast;

my ($protFas,$geneFas,$TEbed,$geneGff,$outdir)=@ARGV;


##1. parse the TAIR blastp result
my $tairPep="/srv/netscratch/dep_mercier/grp_schneeberger/projects/AMPRILdenovo/tair10/TAIR10_pep_20110103_representative_gene_model_updated";
my $tairBlastp="$outdir/evm.ann.prot.blastp.tair10.out";
system("blastp -query $protFas -db $tairPep -evalue 1e-10 -out $tairBlastp -num_threads 10");
parse_blast($tairBlastp,$outdir,"$outdir/tairProt");

##2. parse the TAIR TE related gene blastn result
my $tairTE="/srv/netscratch/dep_mercier/grp_schneeberger/projects/AMPRILdenovo/tair10/TAIR10.TErelated.genes.fasta";
my $tairTEblastn="$outdir/evm.ann.gene.blastn.tair10-TErelated.out";
system("blastn -query $geneFas -db $tairTE -out $tairTEblastn -evalue 1e-5 -num_threads 10");
parse_blast($tairTEblastn,$outdir,"$outdir/tairTEgene");

##3. parse the TE-related protein blastp result
my $TEpeplib="/srv/netscratch/dep_mercier/grp_schneeberger/bin/RepeatMasker/Libraries/RepeatPeps.lib"; 
my $TEblastp="$outdir/evm.ann.prot.blastp.RepeatPebs.out";
system("blastp -query $protFas -db $TEpeplib -evalue 1e-5 -out $TEblastp -num_threads 10");
parse_blast($TEblastp,$outdir,"$outdir/TEprot");
    
##4. intersect TE annotation and gene annotation

open IN,"grep -v '#' $geneGff |grep exon |cut -f 1,4,5,9 |";
open OUT,">$outdir/genes.bed";
while (<IN>) {
  chomp;
  s/model/TU/g;
  my @t=(split /\t/);
  my @tt=(split /\;/,$t[3]);
  my @out=(@t[0..2],@tt);
  my $line=join("\t",@out);
  print OUT "$line\n";
}
close OUT;close IN;

my $overlap = intTEgene("$outdir/genes.bed",$TEbed,"$outdir/gene.TE.intersect.out");
my %overlap = %{$overlap};

#exit;       
## merge results from 1-4 to identity TE-related genes
    
#EVM-Gene	tair10hsp	tairTEhsp	TErelatedHsp	repOverlap%	
open IN,"$outdir/tairProt.blast.besthit.out";
my %tairProt;
while (<IN>) {
  chomp;
  my @t=split /\t/;
  $t[0]=~s/model/TU/g;
  if ($t[1]) {
    $tairProt{$t[0]}=[@t[1,8,9]];	
  }else {
  	$tairProt{$t[0]}=[("None",0,0)];
  }
}    
close IN;

open IN,"$outdir/tairTEgene.blast.besthit.out";
my %tairTE;
while (<IN>) {
  chomp;
  my @t=split /\t/;
  $t[0]=~s/model/TU/g;
  if ($t[1]) {
    $tairTE{$t[0]}=[@t[1,8,9]];	
  }else {
  	$tairTE{$t[0]}=[("None",0,0)];  	
  }
  
}    
close IN;
    
open IN,"$outdir/TEprot.blast.besthit.out";
my %TEprot;
while (<IN>) {
  chomp;
  my @t=split /\t/;
  $t[0]=~s/model/TU/g;
  if ($t[1]) {
    $TEprot{$t[0]} = [@t[1,8,9]];
  }else {
  	$TEprot{$t[0]} = [("None",0,0)];  	
  }
}    
close IN;

open IN,$geneGff;
open OUT1,">$outdir/TE.related.gene.check.out";
open OUT2,">$outdir/evm.all.gene.transposons.gff3";
my %type;
while (<IN>) {
  chomp;
#  next if (!/evm/);  # edit by mg
  next if (!/ID/); 
  my @t = (split /\t/);
  if ($t[2] eq "gene") {
  	my $id=(split /\;/,$t[8])[0];
  	$id=(split /\=/,$id)[1];
  	#print "$id\n";
  	#if ($tair)
  	my @out = (@t[0,3,4,6],$id,@{$tairProt{$id}},@{$tairTE{$id}},@{$TEprot{$id}},$overlap{$id});
  	my $line=join("\t",@out);
  	my $type="U";
  	if ( ($tairProt{$id}[1] > 70) and ($tairProt{$id}[2] > 70) ) {
  	   $type="protein-coding";  		
  	}
  	elsif ( ($tairTE{$id}[1] > 70) and ($tairTE{$id}[2] > 70) ) {
  		$type="TE-related";  			  		
  	}
  	elsif ( ($TEprot{$id}[1] > 70) and ($TEprot{$id}[2] > 70) ) {
  		$type="TE-related";  			  		
  	}
  	elsif ( ($overlap{$id} < 30) and  ($tairProt{$id}[1] > 50) and ($tairProt{$id}[2] > 60)) {
  		$type="protein-coding";  			  		
  	}
  	elsif ( ($overlap{$id} > 30) and  ($tairProt{$id}[1] < 30) and ($tairProt{$id}[2] < 60)) {
  		$type="TE-related";  			  		
  	}
  	elsif ( ($overlap{$id} < 20) and ($tairTE{$id}[1] < 20) and ($tairTE{$id}[2] < 20) and ($TEprot{$id}[1] < 20) and ($TEprot{$id}[2] < 20) ) {
  		$type="protein-coding";  			  		
  	}
  	elsif ( ($overlap{$id} > 30) and ($tairTE{$id}[1] > 30) and ($tairTE{$id}[2] > 30)  ) {
  		$type="TE-related";  			  		
  	}
  	elsif ( ($overlap{$id} > 30) and ($TEprot{$id}[1] > 30) and ($TEprot{$id}[2] > 30)  ) {
  		$type="TE-related";  			  		
  	}
  	elsif ( ($overlap{$id} < 20) and ($tairTE{$id}[1] < 10) and ($tairTE{$id}[2] < 10)  ) {
  		$type="protein-coding";  			  		
  	}
  	elsif ( ($overlap{$id} > 80) and ($tairProt{$id}[1] < 50 ) and ($tairProt{$id}[2] <50 )  ) {
  		$type="TE-related";  			  		
  	}
  	elsif ( ($overlap{$id} > 20) and ($TEprot{$id}[1] > 30) and ($TEprot{$id}[2] > 30) and  ($tairProt{$id}[1] < 50) ) {
  		$type="TE-related";  			  		
  	}
  	elsif (  ($tairProt{$id}[1] > 20 ) and ($tairProt{$id}[2] > 95 )) {
  		$type="protein-coding";
  	}
  	elsif ( $overlap{$id} > 80 ) {
  		$type="TE-related";  			  		
  	}
  	elsif (  ($tairProt{$id}[1] > 20 ) and ($tairProt{$id}[2] > 85 )) {
  		$type="protein-coding";
  	}
  	elsif ( ($tairTE{$id}[1] > 20) and ($tairTE{$id}[2] > 85)  ) {
  		$type="TE-related";  			  		
  	}
  	elsif ($overlap{$id}<30) {
  	  $type="protein-coding";	  		
  	}else {
  		$type="TE-related";
  	}
  	
  	print OUT1 "$line\t$type\n";
  	$type{$id}=$type;
  } 
}
close OUT1;close OUT2;
close IN;    

open IN,$geneGff;
open OUT1,">$outdir/annotation.genes.gff";
open OUT2,">$outdir/annotation.genes.TE.gff";
my $type="";
while (<IN>) {
  chomp;
  next if ! /\d/;
  my @t=(split /\t/);
  if ($t[2] eq "gene" ) {
  	my $id=(split /\;/,$t[8])[0];
  	$id=(split /\=/,$id)[1];
  	if ($type{$id} eq "protein-coding") {
  	  print OUT1 "$_\;Note=protein_coding_gene\n";
  	  print OUT2 "$_\;Note=protein_coding_gene\n";
  	  $type="prot";
  	}else {
  	  print OUT2 "$_\;Note=transposable_element_gene\n";
  	  $type="TE";	
  	}	
  }else {
  	if ($type eq "prot") {
  	  print OUT1 "$_\n";
  	}else {
  	  next if ($t[2] eq "CDS" ) ;
  	}  	
  	print OUT2 "$_\n";
  }  
}  
close IN;close OUT1;close OUT2;  
    

sub intTEgene {
	##intersect the gene annotation gff and TE annotation gff
	my ($gene,$repeat,$out)=@_;
   
	#chr1	5121	5180	ID=Am_100010.1.exon_1	Parent=Am_100010.1
	#chr1	4757	4980	ID=Am_100010.1.exon_2	Parent=Am_100010.1

	#chr1	65941	66159	rnd-5_family-437;DNA/TcMar-Mariner
	#chr1	68545	68671	rnd-4_family-845;LINE/L1
	
	#grep -v '#' chr.all.v1.0.repeats.ann.gff3 |egrep -v 'Low|Simple|Satellite|RNA|Other' |cut -f 1,4,5,9  >chr.all.v1.0.TE.bed
	
	open IN,$gene;
	my %gene;
	while (<IN>) {
	  chomp;
	  my @t=split /\t/;
	  my $id=(split /\=/,$t[4])[1];
	  $gene{$id}+=$t[2]-$t[1]+1;
	}
	close OUT;		
	
	my %overlap;
	open IN,"intersectBed -a $gene -b $repeat -wao|sort -k1,1 -k5,5 -k2,2n -k7,7n |";
	open OUT,">$out";
	my ($g,$overlap,$scaf,$start,$end,$inf,$exonID,$ex);
	$overlap=0;
	my ($ts,$te);
	while (<IN>) {
	  chomp;
	  my @t=split /\t/;
	    
	  my $id=(split /\=/,$t[4])[1]; ##gene id
	  #$id=(split /\./,$id)[1];
	  #pop @tmp;
	  #$id=join("\.",@tmp);
	  
	  ##exon
	  my $exon=(split /\=/,$t[3])[1];
	  
	  if (!$g) {
	    $g=$id;
	    $ex=$exon;
	    $overlap=$t[9]+1 if ($t[9]>0);
	    ($scaf,$start,$end,$inf)=@t[0,1,2,3];  
	    ($ts,$te)=@t[6,7];  
	    
	  }elsif ($g ne $id) {
	    $overlap=sprintf("%.1f",100*$overlap/$gene{$g});
	    print OUT "$scaf\t$g\t$overlap\n";	
	    $overlap{$g}=$overlap;
	    $overlap=0;
	    
	    $g=$id;
	    $ex=$exon;
	    $overlap=$t[9]+1 if ($t[9]>0);
	    ($scaf,$start,$end,$inf)=@t[0,1,2,3];    
	    ($ts,$te)=@t[6,7];
	    
	  }elsif ($ex ne $exon) {
	  	$overlap+=$t[9]+1 if ($t[9]>0);
	  	$ex=$exon;
	  	($ts,$te)=@t[6,7];
	  }else {
	    if ( ($t[6]>=$ts) and ($t[7]<=$te)) {    	
	      next;
	    }elsif ( ($t[6]>=$ts) and ($t[6]<=$te) )  {
	  	   if ($te+1>$t[2]) {
	         next;
	       }else {
	         if ($t[2]<=$t[7]) {
	           $overlap+=$t[2]-($te+1)+1;
	         }else {
	           $overlap+=$t[7]-($te+1)+1;
	         }
	       }
	       $te=$t[7];
	  	}else {
	  	   $overlap+=$t[9]+1;	
	       ($ts,$te)=@t[6,7];
	  	}
	  }
	}
	$overlap=sprintf("%.1f",100*$overlap/$gene{$g});
	$overlap{$g}=$overlap;
	print OUT "$scaf\t$g\t$overlap\n";	
	close OUT;close IN;
	return \%overlap;
}

sub parse_blast {
  my ($blast_result,$outdir,$o)=@_;
  my $in = new Bio::SearchIO(-format => 'blast', -file => "$blast_result",-best_hit_only=>"TRUE");
  open OUT1,">$outdir/$o.blast.besthit.out";
  open OUT2,">$outdir/$o.blast.besthit.mis.gap.out";
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
