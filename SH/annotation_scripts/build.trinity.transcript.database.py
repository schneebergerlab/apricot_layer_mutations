#!/usr/bin/env python
# encoding: utf-8
'''
Created on Jul 26, 2017

@author: jiao
'''
import sys
import os
import getopt
import re
import glob
import subprocess
import time


def main(argv):
    cfgFile = ""
    try:
        opts, args = getopt.getopt(argv,"g:r:b:o:",["refFa=","rnaDir=","bamDir=","outDir=",])  
    except getopt.GetoptError:
        print('build.trinity.transcript.database.py -f <cfg>  ')
        sys.exit(2)
    if len(opts) == 0 :
        print('build.trinity.transcript.database.py -f <cfg>  ')
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print('build.trinity.transcript.database.py -f <cfg>  ')
            sys.exit()
        elif opt in ("-f", "--cfg"):
            cfgFile = arg

    
    #step 1: de novo assembly using Trinity
    
    
    #step 2: genome-guided assembly using Trinity
    

  


if __name__ == "__main__":
   main(sys.argv[1:])


