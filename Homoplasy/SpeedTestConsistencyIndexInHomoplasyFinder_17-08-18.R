#!/usr/bin/env Rscript

#### Load packages ####
suppressMessages(library(homoplasyFinder))

#### Get the input tree and alignment ####

# Get the command line arguments
args = commandArgs(trailingOnly=TRUE)

# Store the arguments
path <- args[1]
treeFile <- paste0(path,args[2])
fastaFile <- paste0(path,args[3])
#path <- "/home/josephcrispell/Desktop/Research/Homoplasy/"
#treeFile <- paste0(path, "example-AFTER_10-08-18.tree")
#fastaFile <- paste0(path, "example_10-08-18.fasta")
createFasta <- TRUE
if(is.na(args[4]) == FALSE){
  createFasta <- as.logical(args[4])
}
createReport <- TRUE
if(is.na(args[5]) == FALSE){
  createReport <- as.logical(args[5])
}
createAnnotatedNewickTree <- TRUE
if(is.na(args[6]) == FALSE){
  createAnnotatedNewickTree <- as.logical(args[6])
}

#### Calculate the consistency index of each site ####
result <- runHomoplasyFinderInJava(treeFile, fastaFile, path, createFasta,
                                   createReport, createAnnotatedNewickTree, verbose=FALSE)
