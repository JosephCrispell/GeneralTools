#!/usr/bin/env Rscript

#### Load packages ####
suppressMessages(library(ape))
suppressMessages(library(phangorn))

#### Get the input tree and alignment ####

# Get the command line arguments
args = commandArgs(trailingOnly=TRUE)

# Store the arguments
path <- args[1]
treeFile <- args[2]
fastaFile <- args[3]
#path <- "/home/josephcrispell/Desktop/Research/Homoplasy/"
#treeFile <- paste0(path, "example-AFTER_10-08-18.tree")
#fastaFile <- paste0(path, "example_10-08-18.fasta")

#### Calculate the consistency index of each site ####

# Read in the tree file
tree <- read.tree(treeFile)

# Read in the sequence file
sequences <- read.phyDat(fastaFile,format="fas")

# Calculate the consistency index
ci <- CI(tree, sequences, sitewise=TRUE)

# Identify inconsistent sites
inconsistent <- (which(!is.nan(ci) & ci < 1))
print(inconsistent)