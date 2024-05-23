#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: add.false.spGenes.pl
#
#        USAGE: ./add.false.spGenes.pl  
#
#  DESCRIPTION: add false halptype specific genes due to unannotation in the other haplotype (indicated by the intersection 
#               between the blastn result and original annotation[from augustus, snap , glimmerhmm])
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Wen-Biao Jiao (), 
# ORGANIZATION: Department of Plant Developmental Biology, Max Planck Institute for Plant Breeding Research
#      VERSION: 1.0
#      CREATED: 03/26/20 18:08:34
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use Bio::SeqIO;
use Bio::Seq;

my ($blnBed,$inGff,$inProt, $indir,$outdir,$p) = @ARGV;
system("mkdir -p $outdir");
# step 1: get the false specific genes and their prot seq
open IN, $blnBed;
my %genes;
while (<IN>) {
  chomp;
  my @t = (split /\t/);
  $genes{$t[9]} = 1;
}
my $seqin = new Bio::SeqIO(-format=>"fasta", -file=>"$inProt");
my %spProt ;
open OUT,">$outdir/$p.missing.other-hap.spGenes.prot.fa";
while (my $seqobj = $seqin->next_seq() ) {
   my $id = $seqobj->id();
   $id = (split /\./,$id)[0];
   #print $id;exit;
   if ($genes{$id} ) {
   	  my $seq = $seqobj->seq();
   	  print OUT ">$id\n$seq\n";
   	  $spProt{$id} = $seqobj->length();  
   }
}
close OUT;
print "to be added genes\t", scalar keys %spProt, "\n";

system("makeblastdb -in $outdir/$p.missing.other-hap.spGenes.prot.fa -dbtype prot");

# step 2: find the potential annotation from abinitio prediction and their prot seq
my %ann;
my @ann = ("augustus", "snap", "glimmerhmm");
my %prot;
foreach my $ann (@ann) {
   open IN,"intersectBed -a $blnBed -b $indir/$ann.bed -f 0.8 -wo |";
   while (<IN>) {
     chomp;	
     my @t = (split /\t/);
     my $id = (split /;/,$t[-2])[0];
     $id = (split /=/,$id)[1];
     push @{$ann{$t[9]}}, $id;
     $prot{$id} = 1;
  }
  close IN;
}
my $n = scalar keys %ann;
my $m = scalar keys %prot;
print "$n genes have $m potential annotation (overlap>0.8) \n";


open OUT,">$outdir/$p.spGenes.the-other.potential.ann.prot.fa";
foreach my $ann (@ann) {
   my $seqin = new Bio::SeqIO(-format=>"fasta", -file=>"$indir/$ann.prot.fa");
   while (my $seqobj = $seqin->next_seq() ) {
       my $id = $seqobj->id();
       if ( $prot{$id}) {
   	     my $seq = $seqobj->seq();
   	     my $ss = substr($seq,0,length($seq)-5);
   	     next if ($ss=~m/\*/);
   	     print OUT ">$id\n$seq\n";  
   	     $prot{$id} = $seqobj->length();
       }
   }
}
close OUT;

# step 3: blast the potential annoation to the other-halptype specific prot seq, find the best hit from the potential annoation
my $blpOut = "$outdir/$p.spGenes.blastn.the-other.prot.out.m6";
my $spProt = "$outdir/$p.spGenes.the-other.potential.ann.prot.fa";
system("blastp -db $outdir/$p.missing.other-hap.spGenes.prot.fa -query $spProt -outfmt 6 -out $blpOut -num_threads 20 ");

open IN, "sort -k2,2 -k1,1 -k9,9n $blpOut |";
my %best; my %add;
while (<IN>) {
   chomp;
   my @t = (split /\t/);
   if ( $best{$t[1]}) {
   	  if ($t[-1] > $best{$t[1]}[4]) {
   	  	 my $iden = $t[2];
   	     my $cov1 = $t[3]/$prot{$t[0]};
   	     my $cov2 = $t[3]/$spProt{$t[1]};
   	  	 @{$best{$t[1]}} = ($t[0], $iden,$cov1,$cov2, $t[-1] ); 
   	  }
   }else {
   	  my $iden = $t[2];
   	  my $cov1 = $t[3]/$prot{$t[0]};
   	  my $cov2 = $t[3]/$spProt{$t[1]};
   	  if ($iden > 80 and $cov1 > 0.80 and $cov2 > 0.80 ) {
   	  	@{$best{$t[1]}} = ($t[0], $iden,$cov1,$cov2, $t[-1] ); 	
   	  }
      
   }	
}
close IN;
foreach my $g (sort keys %best) {
   $add{$best{$g}[0]} = 1;
}
print "## can be added gene model\t", scalar keys %best, "\n";



# step 4: add gene models
system("cat $inGff >$outdir/$p.add.missing.ann.gff");
open OUT,">>$outdir/$p.add.missing.ann.gff";
my %lines;
foreach my $ann (@ann) {
   open IN,"$indir/$ann.gff";
   my $flag = 0; my $id ;
   while (<IN>) {
     chomp;	
     my @t = (split /\t/);
     if ($t[2] eq "gene") {
     	$id = (split /;/,$t[8])[0];
        $id = (split /=/,$id)[1];
        if ($add{$id} ) {
           	push @{$lines{$id}}, $_; 
           	$flag = 1;
           	print OUT "$_\n";
        }else {
           $flag = 0;	
           #$id = "";
        }
     }else {
     	if ($flag == 1) {
     	   push @{$lines{$id}}, $_;	
     	   print OUT "$_\n";
     	}
     }
  }
  close IN;
}
close OUT;

#step 5:  sort the gene model 
open IN,"$outdir/$p.add.missing.ann.gff";
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
    #print keys %genes2;
  }
}
$genes1{$ch}{$id}=$start;
$genes2{$ch}{$id}=[@tmp];

#print keys %genes2;
#print $ch,"\n";

my $out = "$outdir/$p.add.missing.ann.sorted.gff";
open OUT,">$out";
my %newGenes;
foreach my $chr (sort keys %genes1) {
  foreach my $id (sort {$genes1{$chr}{$a}<=>$genes1{$chr}{$b}} keys %{$genes1{$chr}}) {
    if (!$genes2{$chr}{$id}) {
      print  "$chr\t$id\n";
      exit;
    }
    my $line=join("\n",@{$genes2{$chr}{$id}});
    print OUT "$line\n";
    push @{$newGenes{$chr}}, $id;
  }
}
close OUT;

#step 6: replace the added model with new gene ID
my $oo = $p;
$p = uc($p);
my %replID ;

foreach my $chr (sort keys %newGenes) {
  my @tmp = @{$newGenes{$chr}};
  foreach my $i (0..$#tmp) {
     if (substr($newGenes{$chr}[$i],0,4) ne "$p-") {
       	my $j = $i-1;
       	my $pre = -1;
       	while ($j>=0) {
       	  my $preID = $newGenes{$chr}[$j];
       	  if (substr($preID, 0, 4) eq "$p-") {
       	  	$pre = $j;
       	    last;	
       	  }else {
       	  	 $j--;
       	  }
       	}
       $j = $i + 1;
       my $nex = $#tmp + 1;
       while ($j<=$#tmp) {
       	  my $nexID = $newGenes{$chr}[$j];
       	  if (substr($nexID, 0, 4) eq "$p-") {
       	  	$nex = $j;
       	    last;	
       	  }else {
       	  	 $j++;
       	  }
       }
       print "index: $pre\t$nex\n";
       my ($pid,$nid);
       if ($pre == -1) {
       	  $pid = 10000;
       	  $nid = 10010;
       	  print "to be added at the start  $pid\t$nid\t$newGenes{$chr}[$i]\n";
       }
       elsif ($nex==$#tmp + 1) {
         $pid =  substr($newGenes{$chr}[$pre],6,length($newGenes{$chr}[$pre])-6);
         $nid = $pid + 20;
         $nex = $i + 1;
         print "to be added at the end  $nid\t$pid\t$newGenes{$chr}[$i]\n";
       }else {
       	 $pid = substr($newGenes{$chr}[$pre],6,length($newGenes{$chr}[$pre])-6);
       	 $nid =  substr($newGenes{$chr}[$nex],6,length($newGenes{$chr}[$nex])-6);
       	 print "id  $pid\t$nid\t$newGenes{$chr}[$i]\n";
       }
       
       my $flag =0; my $span = 5;
       
       if (abs($nid - $pid) <= $nex-$i ) {
       	   my $kk = $nex - $i ;
       	   print "$kk more genes to added between $chr $nid and $pid\n";	
       	   $span = int(10*($nid-$pid)/($nex-$i+1))*10;	
       	   $flag = 1;
       }else {
       	   $span = int(abs($nid-$pid)/($nex-$i+1));	
       }
       $span = 5 if ($span == 10);
       print "span $span\n";
       foreach my $k ($i..$nex-1) {
         	if (!$flag) {
         	  my $newID = $pid + $span*($k-$i+1);
         	  $replID{$tmp[$k]} = $p ."-". substr($chr,3,1) ."G". $newID;
         	  $newGenes{$chr}[$k] = $p ."-". substr($chr,3,1) ."G". $newID;
         	}else {
         	  my $newID = $pid + $span*($k-$i+1);
         	  $newID = $newID . "5";
         	  $replID{$tmp[$k]} = $p ."-". substr($chr,3,1) ."G". $newID;
         	  $newGenes{$chr}[$k] = $p ."-". substr($chr,3,1) ."G". $newID;
         	} 	
       }	      
     }	
  }
}

$out = "$outdir/$oo.add.missing.ann.sorted.gff";
open IN,"$out";
my $out2 = "$outdir/$oo.add.missing.ann.sorted.newID.gff";
open OUT,">$out2";
open OUT2,">$outdir/$oo.add.newID.translation.txt";
open OUT3,">$outdir/$oo.new.protein-coding.genes.gff";
my $k = 0; my $j = 0;
my $flag = 0; my $isTE = 0;
my $id = "";
while (<IN>) {
  chomp;
  my @t=split /\t/;
  if ($t[2] eq "gene") {
     $id = (split /;/,$t[8])[0];
     $id = (split /=/,$id)[1];
     if ( $replID{$id}) {
       my @tmp = (split /;/,$t[8]);
       shift @tmp;
       my $tmp = join(";",@tmp);
       $t[8] = "ID=$replID{$id}\;Note=protein_coding_gene";
       my $line = join("\t",@t);
       print OUT "$line\n";	
       print OUT2 "$id\t$replID{$id}\n";
       print OUT3 "$line\n";	
       $k = 1; $j = 1;
       $flag = 1;
     }else {
     	$flag = 0;
     	print OUT "$_\n";
     	if ($t[8]=~/trans/) {
     	  $isTE = 1;	
     	}else {
     	   $isTE = 0;
     	   print OUT3 "$_\n";		
     	}
     }
  }else {
     if ($flag == 0) {
       print OUT "$_\n";	
       if ($isTE == 0) {
         print OUT3 "$_\n";		
       }
     }else {
     	if ($t[2] eq "mRNA") {
     	   $t[8] = "ID=$replID{$id}\.1\;Parent=$replID{$id}";
     	   my $line = join("\t",@t);
           print OUT "$line\n";	
           print OUT3 "$line\n";	
     	}elsif ($t[2] eq "exon") {
     	   $t[8] = "ID=$replID{$id}\.1.exon.$k\;Parent=$replID{$id}\.1";
     	   my $line = join("\t",@t);
           print OUT "$line\n";
           print OUT3 "$line\n";	
           $k ++;
     	}elsif ($t[2] eq "CDS") {
     		$t[8] = "ID=$replID{$id}\.1.cds.$j\;Parent=$replID{$id}\.1";
     	   my $line = join("\t",@t);
           print OUT "$line\n";
           print OUT3 "$line\n";	
           $j ++ ;
     	}
     }
  }
}
close OUT;close OUT2;close OUT3;

open OUT,">$outdir/ungr.gene.can.be.added.txt";
foreach my $g (sort keys %best) {
	#ungr added newID iden cov cov2
   print OUT "$g\t$best{$g}[0]\t$best{$g}[1]\t$best{$g}[2]\t$best{$g}[3]\t$replID{$best{$g}[0]}\n";
}
close OUT;

