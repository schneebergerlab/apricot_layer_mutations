#TE annotation

1) RepeatModeler
#Source Distribution Installation

## Prerequisites

#Perl Available at http://www.perl.org/get.html. Developed and tested with version 5.8.8.
/srv/netscratch/dep_coupland/grp_schneeberger/bin/perl/perl-5.26.1/local/bin/perl

#RepeatMasker & Libraries Developed and tested with open-4.0.9. The program is available at http://www.repeatmasker.org/RMDownload.html and is distributed with Dfam - an open database of transposable element families.
/srv/netscratch/dep_coupland/grp_schneeberger/bin/RepeatModeler/RepeatMasker/RepeatMasker

#RECON - De Novo Repeat Finder, Bao Z. and Eddy S.R. Developed and tested with our patched version of RECON ( 1.08 ). The 1.08 version fixes problems with running RECON on 64 bit machines and supplies a workaround to a division by zero bug along with some buffer overrun fixes. The program is available at: http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz. The original version is available at http://eddylab.org/software/recon/.
/srv/netscratch/dep_mercier/grp_schneeberger/bin/RECON/RECON-1.08/bin/

#RepeatScout - De Novo Repeat Finder, Price A.L., Jones N.C. and Pevzner P.A. Developed and tested with our multiple sequence version of RepeatScout ( 1.0.6 ). This version is available at http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
/home/jiao/bin/RepeatScout

#TRF - Tandem Repeat Finder, G. Benson et al. You can obtain a free copy at http://tandem.bu.edu/trf/trf.html RepeatModeler requires version 4.0.9 or higher.
/home/jiao/bin/trf

#RMBlast - A modified version of NCBI Blast for use with RepeatMasker and RepeatModeler. Precompiled binaries and source can be found at http://www.repeatmasker.org/RMBlast.html We recommend using 2.9.0-p2 or higher.
/srv/netscratch/dep_mercier/grp_schneeberger/bin/RepeatModeler/rmblast-2.10.0$ cp bin/* /home/jiao/bin/


#Optional. Additional search engine:
#ABBlast - Sequence Search Engine, W. Gish et al. Developed and tested with 2.0 [04-May-2006]. The latest versions of ABBlast may be downloaded from: http://blast.advbiocomp.com/licensing/


#Optional. Required for running LTR structural search pipeline:
#LtrHarvest - The LtrHarvest program is part of the GenomeTools suite. We have developed this release of RepeatModeler on GenomeTools version 1.5.9 available for download from here: http://genometools.org/pub/ NOTE: use the "make threads=yes" build options to enable multi-threaded runs.
/srv/netscratch/dep_mercier/grp_schneeberger/bin/genometools/genometools-1.5.9/bin/bin/gt

#Ltr_retriever - A LTR discovery post-processing and filtering tool. We recommend using version 2.6 or higher from here: https://github.com/oushujun/LTR_retriever/releases
/srv/netscratch/dep_mercier/grp_schneeberger/bin/LTR_retriever/LTR_retriever-2.8/bin
/srv/netscratch/dep_mercier/grp_schneeberger/bin/cdhit/cd-hit-v4.8.1-2019-0228 #cdhit
/srv/netscratch/dep_coupland/grp_schneeberger/bin/hmmer/hmmer-3.2.1/local/bin #hmmer


#MAFFT - A multiple sequence alignment program. We developed and tested RepeatModeler using mafft version 7.407. Please use this verison or higher from here: https://mafft.cbrc.jp/alignment/software/
/opt/share/software/bin/

#CD-HIT - A sequence clustering package. We developed and tested RepeatModeler using version 4.8.1. Please use this version or higher from: http://weizhongli-lab.org/cd-hit/
/srv/netscratch/dep_mercier/grp_schneeberger/bin/cdhit/cd-hit-v4.8.1-2019-0228

#Ninja - A tool for large-scale neighbor-joining phylogeny inference and clustering. We developed and tested RepeatModeler using Ninja version "0.95-cluster_only". Please obtain a copy from: https://github.com/TravisWheelerLab/NINJA/releases/tag/0.95-cluster_only
/srv/netscratch/dep_mercier/grp_schneeberger/bin/RepeatModeler/NINJA-0.95-cluster_only/NINJA


Installation

Obtain the source distribution
From repeatmasker.org:
Latest Version Released 1/9/2019: RepeatModeler-2.0.1.tar.gz
Previous Version Released 11/15/2019: RepeatModeler-2.0.tar.gz
Uncompress and expand the distribution archive:

Typically:
#tar -zxvf RepeatModeler-open-#.#.#.tar.gz
or
#gunzip RepeatModeler-open-#.#.#.tar.gz
#tar -xvf RepeatModeler-open-#.#.#.tar



Automatic:
#Run the "configure" script interactively with prompts for each setting:

   perl ./configure
#Run the "configure" script with supplied paramters:

   perl ./configure -rscout_dir .. -recon_dir ..
By Hand:
Edit the configuration file "RepModelConfig.pm"
Dynamically:
Use the "configuration overrides" command line options with the RepeatModeler programs. e.g:

   ./RepeatModeler -rscout_dir .. -recon_dir .. 

http://www.repeatmasker.org/RepeatModeler/

################# 
#EDTA
https://github.com/oushujun/EDTA

Step by step installation using conda

conda create -n EDTA
conda activate EDTA
conda config --env --add channels anaconda --add channels conda-forge --add channels bioconda
conda install -n EDTA -y cd-hit repeatmodeler muscle mdust blast java-jdk perl perl-text-soundex multiprocess regex tensorflow=1.14.0 keras=2.2.4 scikit-learn=0.19.0 biopython pandas glob2 python=3.6 tesorter genericrepeatfinder genometools-genometools ltr_retriever ltr_finder
git clone https://github.com/oushujun/EDTA
./EDTA/EDTA.pl


###### RNAseq #########
# Trinity
/srv/netscratch/dep_mercier/grp_schneeberger/bin/Trinity/trinityrnaseq-v2.9.1/Trinity #run under base ENV

# salmon


#####Maker######
(EDTA) jiao@build-stretch:/srv/netscratch/dep_mercier/grp_schneeberger/bin/maker/maker-3.01 #install under EDTA env

###PASA
# https://sr-c.github.io/2018/07/15/pasa-install/




