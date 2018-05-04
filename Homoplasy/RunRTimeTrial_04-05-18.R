library(devtools)
install_github("JosephCrispell/homoplasyFinder")
library(homoplasyFinder)
library(ape)


system.time(run(100))

run <- function(nSequences){
  path <- "/home/josephcrispell/Desktop/Research/Homoplasy/TimeTrial/"
  fastaFile <- paste(path, "Example_", 2* nSequences, "-", nSequences, "_27-04-18.fasta", sep="")
  treeFile <- paste(path, "Example_", 2* nSequences, "-", nSequences, "_27-04-18.tree", sep="")
  
  sequences <- read.dna(fastaFile, skip=1, format="fasta")
  tree <- read.tree(treeFile)
  
  homoplasyFinder(tree, sequences, verbose=FALSE)
}
