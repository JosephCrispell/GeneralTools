#### Preparation ####

# Load libraries
library(phyloHelpeR) # Build RAXML tree and tangle plot devtools::install_github("JosephCrispell/phyloHelpeR")
library(phytools) # Rotate node on phlogeny
library(ape) # ladderize tip - reordering tips to match
library(ade4) # mantel test to compare matrices

# Set the path
path <- file.path("~", "storage", "Research", "RepublicOfIreland", "Mbovis", "Wicklow")
  
# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Build the phylogeny for the fasta with PE/PPE regions ####

# Build the phylogeny using the FASTA with PE/PPE regions
fastaFilePEPPE <- file.path(path, "vcfFiles", "sequences_Prox-10_29-03-2020.fasta")
treeWithPEPPE <- runRAXML(fastaFilePEPPE, date="29-03-20", path, nThreads=10, outgroup="\\>Ref-1997")

# Read in the phylogeny built without the PE/PPE regions
fastaFile <- file.path(path, "vcfFiles", "sequences_Prox-10_29-03-2020.fasta")
tree <- runRAXML(fastaFile, date="21-10-19", path, nThreads=10, outgroup="\\>Ref-1997")

# Remove NI isolates and Reference
treeWithPEPPE <- drop.tip(treeWithPEPPE, treeWithPEPPE$tip.label[grepl(treeWithPEPPE$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])

# Parse the tip labels
treeWithPEPPE$tip.label <- editTipLabels(treeWithPEPPE$tip.label)
tree$tip.label <- editTipLabels(tree$tip.label)

# Re-scale the branch lengths to be SNPs
treeWithPEPPE$edge.length <- treeWithPEPPE$edge.length * getNSitesInFASTA(fastaFilePEPPE)
tree$edge.length <- tree$edge.length * getNSitesInFASTA(fastaFilePEPPE)

# Re-order the trees so that they match as much as possible
treeWithPEPPE <- ladderize(treeWithPEPPE, right=FALSE)
tree <- ladderize(tree, right=FALSE)

#### Plot a tangle plot to illustrate the differences between the phylogenies ####

tanglePlot(tree, treeWithPEPPE, connectingLine.col="red", connectingLine.lwd=1.5, connectingLine.lty=2, offsetProp=0.065)

phyloHelpeR::addSNPScale(position="topright")

#### Compare the patristic distances ####

# Read in the fasta file and generate genetic distance matrices
sequences <- read.dna(fastaFilePEPPE, format="fasta")
distancesWithPEPPE <- dist.dna(sequences)
sequences <- read.dna(fastaFile, format="fasta")
distances <- dist.dna(sequences)

# Calculate the correlation between the matrices (check https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/)
results <- mantel.rtest(distancesWithPEPPE, distances, nrepet = 9999)

plot(distances, distancesWithPEPPE, pch=19, col=rgb(0,0,0, 0.1), cex=0.5)

#### FUNCTIONS ####

getNSitesInFASTA <- function(fastaFile){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get first line of file
  firstLine <- readLines(connection, n=1)
  
  # Close file connection
  close(connection)
  
  # Get the number of sites used in the FASTA file from first line
  nSites <- as.numeric(strsplit(firstLine, " ")[[1]][2])
  
  return(nSites)
}

editTipLabels <- function(tipLabels){
  
  # Initialise a vector to store the new labels
  output <- c()
  
  # Examine each tip label
  for(label in tipLabels){
    
    # Check if starts with ">"
    if(grepl(label, pattern="^>")){
      
      # Remove the ">" prefix
      label <- substr(label, 2, nchar(label))
    }
    
    # Split the label and retain first part
    label <- strsplit(label, split="_")[[1]][1]
    
    # Store the new label
    output[length(output) + 1] <- label
  }
  
  return(output)
}
