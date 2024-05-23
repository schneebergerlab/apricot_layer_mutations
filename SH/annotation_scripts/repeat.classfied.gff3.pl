#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: repeat.classfied.gff3.pl
#
#        USAGE: ./repeat.classfied.gff3.pl  
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
#      CREATED: 05/22/2014 03:48:03 PM
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;
use File::Basename;

my ($gff,$rmout,$out,$stat)=@ARGV;

open IN,$rmout;
my %repeat;
while (<IN>) {
  chomp;
  next if (!/\(\d/);
  my $tmp=" ".$_;
  my @line=split /\s+/,$tmp;
  #print "$line[10]\t$line[11]\n";exit;
  $repeat{$line[5]}{$line[6]}="$line[10]\t$line[11]";
}

open IN,$gff;
my $n=0;
open OUT,">$out";
open OUT2,">$stat";
my %class=();
while (<IN>) {
  chomp;
  if (/#/) {
  	print OUT "$_\n";
  	next;
  }
  $n++;
  my @line=split /\t/;
  my $infor=$repeat{$line[0]}{$line[3]};
  my ($fa,$class)=split /\t/,$infor;
  if ($class=~m/\//) {
  	my ($c,$l)=split /\//,$class;
  	$infor="ID=".$line[0].".TE.".$n."\;Subfamily=$fa\;Family=$l\;Class=$c";
  	$class{$c}{$l}{"num"}++;
  	my $leng=$line[4]-$line[3]+1;
  	$class{$c}{$l}{"leng"}+=$leng; 	
  } 
  else {
  	$infor="ID=".$line[0].".TE.".$n."\;Subfamily=$fa\;Family=$class\;Class=$class";
  	$class{$class}{$class}{"num"}++;
  	my $leng=$line[4]-$line[3]+1;
  	$class{$class}{$class}{"leng"}+=$leng;
  }
  my $line=join("\t",@line[0..7],$infor);
  print OUT "$line\n";
}
my $totalL;my $totalN;
foreach my $class (sort keys %class) {
  my %tmp=%{$class{$class}};
  my $leng=0;my $num=0;
  foreach my $fa (sort keys %tmp) {
  	print OUT2 "$class\t$fa\t$class{$class}{$fa}{'num'}\t$class{$class}{$fa}{'leng'}\n";
  	$num+=$class{$class}{$fa}{'num'};
  	$leng+=$class{$class}{$fa}{'leng'};
  	$totalL+=$class{$class}{$fa}{'leng'};
  	$totalN+=$class{$class}{$fa}{'num'};
  }
  print OUT2 "$class\t\-\t$num\t$leng\n"; 
}
print OUT2 "Total\t\-\t$totalN\t$totalL\n";
