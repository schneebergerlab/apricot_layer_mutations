#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: noncoding.infernal.output.parser.pl
#
#        USAGE: ./noncoding.infernal.output.parser.pl  
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
#      CREATED: 06/17/2014 01:39:21 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use File::Basename;

my ($in,$outdir)=@ARGV;

##gff3 output format:

open IN,"cat $in | tr -s ' ' |";
#open OUT,">$outdir/noncoding.gene.gff3";
system("mkdir -p $outdir/tmp");
open OUT1,"|sort -k1,1 -k4,4n -k5,5n >$outdir/tmp/tRNA.tmp.gff3";
open OUT2,"|sort -k1,1 -k4,4n -k5,5n >$outdir/tmp/snoRNA.tmp.gff3";
open OUT3,"|sort -k1,1 -k4,4n -k5,5n >$outdir/tmp/mirRNA.tmp.gff3";
open OUT4,"|sort -k1,1 -k4,4n -k5,5n >$outdir/tmp/snRNA.tmp.gff3";
open OUT5,"|sort -k1,1 -k4,4n -k5,5n >$outdir/tmp/rRNA.tmp.gff3";
open OUT6,"|sort -k1,1 -k4,4n -k5,5n >$outdir/tmp/other.tmp.gff3";
my ($n1,$n2,$n3,$n4,$n5,$n6);
while (<IN>) {
  chomp;
  next if (/\#/);
  next if (/\?/);
  my ($name,$acc,$seq,$start,$end,$strand,$trunc,$evalue,$desc)=(split /\s+/)[0,1,2,7,8,9,10,15,17];
  #print "$name,$acc,$seq,$start,$end,$strand,$trunc,$desc\n";exit;
  if ($start>$end) {
  	($start,$end)=($end,$start);
  }
  if ($trunc eq "no") {
  	if (($name=~m/tRNA/) or ($name=~m/tmRNA/)) {  ##
  	  $n1++;
  	  print OUT1 "$seq\tInfernal_v1_1\ttRNA\t$start\t$end\t$evalue\t$strand\t\.\tID=tRNA.$n1\;targetName=$name\;Pfam=$acc\n";	
  	}
  	elsif ($name=~m/sno/i) {
  	  $n2++;
  	  print OUT2 "$seq\tInfernal_v1_1\tsnoRNA\t$start\t$end\t$evalue\t$strand\t\.\tID=snoRNA.$n2\;targetName=$name\;Pfam=$acc\n";
  	}
  	elsif ($name=~m/mir/i) {
  	  $n3++;
  	  print OUT3 "$seq\tInfernal_v1_1\tmirRNA\t$start\t$end\t$evalue\t\.\t\.\tID=mirRNA.$n3\;targetName=$name\;Pfam=$acc\n";	
  	}
  	elsif ($name=~m/U\d/) {
  	  $n4++;
  	  print OUT4 "$seq\tInfernal_v1_1\tsnRNA\t$start\t$end\t$evalue\t$strand\t\.\tID=snRNA.$n4\;targetName=$name\;Pfam=$acc\n";	
  	}
  	elsif ($name=~m/rRNA/) {
  	  $n5++;
  	  print OUT5 "$seq\tInfernal_v1_1\trRNA\t$start\t$end\t$evalue\t$strand\t\.\tID=rRNA.$n5\;targetName=$name\;Pfam=$acc\n";	
  	}
  	else {
      $n6++;
  	 print OUT6 "$seq\tInfernal_v1_1\tothers\t$start\t$end\t$evalue\t$strand\t\.\tID=others.$n6\;targetName=$name\;Pfam=$acc\n";	
  	}  	
  }
}
close OUT1;close OUT2;close OUT3;close OUT4;close OUT5; close OUT6;
my @rna = ("tRNA","snoRNA","mirRNA","snRNA","rRNA","other");
foreach my $rna (@rna) {
  removeOverlap("$outdir/tmp/$rna.tmp.gff3","$outdir/$rna.gff3");	
}

#system("mv $outdir/*.tmp.gff3 $outdir/tmp");
system("rm $outdir/noncoding.gene.gff3") if (-e "$outdir/noncoding.gene.gff3");
system("cat $outdir/*.gff3 |grep chr |sort -k1,1 -k4,4n >  $outdir/tmp/noncoding.gene.tmp.gff3");
system("cat $outdir/*.gff3 |grep -v chr |sort -k1,1 -k4,4n >>  $outdir/tmp/noncoding.gene.tmp.gff3");


open IN,"$outdir/tmp/noncoding.gene.tmp.gff3";
open OUT,">$outdir/noncoding.gene.gff3";
while (<IN>) {
	chomp;
	my @t=split /\t/;
	if (/\d_/) {
	  my $n=(split /\_/,$t[0])[1];
	  $t[3]+=$n*10000000;
	  $t[4]+=$n*10000000;
	  $t[0]=(split /\_/,$t[0])[0];
	}
	my $line=join("\t",@t);
	print OUT "$line\n";
}
close IN;
close OUT;





sub removeOverlap {
	my ($in,$out)=@_;
	#chr1_0  Infernal_v1_1   snRNA   3462002 3462217 .   +   .   ID=snRNA.10;Name=U3;Description=-;Pfam=RF00012
	#chr1_0  Infernal_v1_1   snRNA   3462002 3462218 .   +   .   ID=snRNA.1;Name=Plant_U3;Description=-;Pfam=RF01847
	my ($ctg,$s,$e)=("",0,0);
	open IN,"$in";
	#print "$in\n";exit;
	my @group;
	my @tmp;
	while (<IN>) {
	  chomp;
	  my @t = split /\t/;
	  if (!$ctg) {
	  	($ctg,$s,$e)= @t[0,3,4];
	  	push @tmp,$_;
	  }elsif ($ctg ne $t[0]) {
	  	($ctg,$s,$e)= @t[0,3,4];
	  	push @tmp,$_;
	  }else {
	  	if ($t[3] < $e) {	  		
	  	  push @tmp,$_;
	  	  $e=$t[4] if ($e < $t[4]);
	  	}else {
	  	  push @group,[@tmp];
	  	  @tmp=();
	  	  push @tmp,$_;
	  	  $e=$t[4];	
	  	}
	  }	  	
	}
	push @group,[@tmp];
	close IN;
	
	open OUT,">$out";
	foreach my $g (@group) {
	  my @t = @{$g};
	  if ($#t==0) {
	  	print OUT "$t[0]\n";
	  }else {
	  	my $line=$t[0];
	  	my $eval=(split /\t/,$line)[5];
	  	foreach my $k (1..$#t) {
	  	  my $ev=(split /\t/,$t[$k])[5];
	  	  if ($ev < $eval) {
	  	  	$line=$t[$k];
	  	  }
	  	}
	  	print OUT "$line\n";
	  }
	}
	close OUT;
				
}


