#### Preparation ####

# Load libraries
library(phyloHelpeR) # Build RAXML tree and tangle plot
library(phytools) # Rotate node on phlogeny

# Set the path
path <- file.path("~", "storage", "Research", "RepublicOfIreland", "Mbovis", "Wicklow")
  
# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

#### Build the phylogeny for the fasta with PE/PPE regions ####

# Build the phylogeny using the FASTA with PE/PPE regions
fastaFile <- file.path(path, "vcfFiles", "sequences_Prox-10_29-03-2020.fasta")
treeWithPEPPE <- runRAXML(fastaFile, date="29-03-20", path, nThreads=10, outgroup="\\>Ref-1997")

# Read in the phylogeny built without the PE/PPE regions
fastaFile <- file.path(path, "vcfFiles", "sequences_Prox-10_29-03-2020.fasta")
tree <- runRAXML(fastaFile, date="21-10-19", path, nThreads=10, outgroup="\\>Ref-1997")

# Remove NI isolates and Reference
treeWithPEPPE <- drop.tip(treeWithPEPPE, treeWithPEPPE$tip.label[grepl(treeWithPEPPE$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])

# Edit the tip labels
treeWithPEPPE$tip.label <- editTipLabels(treeWithPEPPE$tip.label)
tree$tip.label <- editTipLabels(tree$tip.label)

# Reorder main clades in tree to help visualisation
treeWithPEPPE <- rotateNodes(treeWithPEPPE, nodes=c(47,56))
tree <- rotateNodes(tree, nodes=c(47,56))


#### Plot a tangle plot to illustrate the differences between the phylogenies ####

tanglePlot(tree, treeWithPEPPE)

#### FUNCTIONS ####

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
