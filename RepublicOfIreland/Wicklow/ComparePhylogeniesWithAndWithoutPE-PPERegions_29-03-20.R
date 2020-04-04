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

# Open a PDF
pdf(file.path(path, "Figures", paste0("tanglePlot_", date, ".pdf")))

# Plot the two phylogenies facing one another and highlight differences
tanglePlot(tree, treeWithPEPPE, connectingLine.col="red", connectingLine.lwd=1.5, connectingLine.lty=2, offsetProp=0.08)

# Add a SNP scale
addSNPScale(position="topright", lineWidth=1.5)

# CLose the PDF
dev.off()

#### Compare the genetics distances ####

# Read in the fasta file and generate genetic distance matrices
distancesWithPEPPE <- dist.dna(read.dna(fastaFilePEPPE, format="fasta"), model="N")
distances <- dist.dna(read.dna(fastaFile, format="fasta"), model="N")

# Calculate the correlation between the matrices (check https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/)
results <- mantel.rtest(distancesWithPEPPE, distances, nrepet = 9999)

## Plot the genetics distances against one another

# Open a PDF
pdf(file.path(path, "Figures", paste0("ComparingGeneticDistances_", date, ".pdf")))

# Turn off scientific notation for lare/small numbers
options(scipen=10) # Default is zero

# Create y = x set of values
range <- range(distances)
values <- seq(range[1], range[2])

# Create the initial plot with X=Y line
plot(values, values, type="l", lty=2, lwd=2, col="red", bty="n", las=1,
     xlab="Number of SNVs WITHOUT PE/PPE regions",
     ylab="Number of SNVs WITH PE/PPE regions",
     main="Comparing genetic distances calculated with/without PE/PPE regions")

# Add the genetic distances
points(distances[lower.tri(distances)], distancesWithPEPPE[lower.tri(distancesWithPEPPE)], pch=19, col=rgb(0,0,0, 0.05), cex=2)

# Add the mantel.test results
legend("topleft", legend=paste0("Mantel test observation: ", results$obs, " (p = ", round(results$pvalue, digits=5), ")"), bty="n")

# Close the PDF
dev.off()

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
