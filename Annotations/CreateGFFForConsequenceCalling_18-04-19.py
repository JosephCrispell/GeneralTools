#!/usr/bin/env python

# Author Damien Farrell
# Date: 09-02-18
# Purpose: bcftools csq takes a gff file but it doesn't work on the NCBI version. This script creates a GFF from a GB file that will work

# Command line structure:
# python CreateGFFForConsequenceCalling_DATE.py annotations.gb output.gff
#   annotations.gb      Annotations in the Genbank (full) format - download from NCBI
#   output.gff          Annotations in a GFF format that will work with bcftools csq

###################
# Import packages #
###################

import sys # Command line arguments
import BCBio # GFF file object
import Bio.SeqFeature # Feature holders
from Bio import SeqIO

#############
# FUNCTIONS #
#############

def GFF_bcftools_format(in_file, out_file):
    """Convert a bacterial genbank file from NCBI to a GFF3 format that can be used in bcftools csq.
    see https://github.com/samtools/bcftools/blob/develop/doc/bcftools.txt#L1066-L1098.
    Args:
        in_file: genbank file
        out_file: name of GFF file
    """

    from BCBio import GFF
    in_handle = open(in_file)
    out_handle = open(out_file, "w")
    from Bio.SeqFeature import SeqFeature
    from Bio.SeqFeature import FeatureLocation
    from copy import copy

    for record in SeqIO.parse(in_handle, "genbank"):
        #make a copy of the record as we will be changing it during the loop
        new = copy(record)
        new.features = []
        #loop over all features
        for i in range(0,len(record.features)):          
            feat = record.features[i]
            q = feat.qualifiers
            #remove some unecessary qualifiers
            for label in ['note','translation','product','experiment']:
                if label in q:
                    del q[label]
            if(feat.type == "CDS"):
                #use the CDS feature to create the new lines
                tag = q['locus_tag'][0]
                q['ID'] = 'CDS:%s' %tag
                q['Parent'] = 'transcript:%s' %tag
                q['biotype'] = 'protein_coding'

                #create mRNA feature
                m = SeqFeature(feat.location,type='mRNA',strand=feat.strand)
                q = m.qualifiers
                q['ID'] = 'transcript:%s' %tag
                q['Parent'] = 'gene:%s' %tag
                q['biotype'] = 'protein_coding'
                new.features.append(m)

            elif(record.features[i].type == "gene"):
                #edit the gene feature
                q=feat.qualifiers
                q['ID'] = 'gene:%s' %q['locus_tag'][0]
                q['biotype'] = 'protein_coding'
                if 'gene' in q:
                    q['Name'] = q['gene']
            new.features.append(feat)
        #write the new features to a GFF                                      
        GFF.write([new], out_handle)
        return

###########################
# Editing annotation file #
###########################

# Get the command line arguments
annotationFile = sys.argv[1]
outputFile = sys.argv[2]

# Edit the annotation file and create new output
GFF_bcftools_format(annotationFile, outputFile)
