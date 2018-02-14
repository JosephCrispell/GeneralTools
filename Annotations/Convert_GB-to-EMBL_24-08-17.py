#!/usr/bin/env python

# Author Damien Farrell
# Date: 09-02-18

import sys,os
from Bio import SeqIO

def convert_sequence_format(infile, outformat='embl'):
    """convert sequence files using SeqIO"""

    informat = os.path.splitext(infile)[1][1:]
    if informat == 'fa':
        informat = 'fasta'
    print ('input format: %s' %informat)
    print ('output format: %s' %outformat)
    outfile = os.path.splitext(infile)[0]+'.'+outformat
    count = SeqIO.convert(infile, informat, outfile, outformat)
    print ("Converted %i records" %count)
    return

infile = sys.argv[1]
convert_sequence_format(infile)