#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: get.TE.gff.pl
#
#        USAGE: ./get.TE.gff.pl  
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
#      CREATED: 03/23/20 22:29:01
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use utf8;


my ($in,$out) =@ARGV;

open IN,$in;
open OUT,">$out";
my $id;
while (<IN>) {
	chomp;
	my @t = (split /\t/);
	if ($t[2] eq "gene") {
		if  ($t[8]=~/trans/ ) {
			$id = (split /;/,$t[8])[0];
			$id = (split /=/,$id)[1];
			#print $id;exit;
			print OUT "$_\n";
		}else {
			$id = "";
		}
	}else {
		if ($id) {
			print OUT "$_\n";
		}else {
			next;
		}
	}
	
}
close IN;close OUT;


