#### PACKAGES ####
library(ape)
library(phangorn)
library(phytools)
library(geiger)


#### READ FASTA ####

# Read in the FASTA file
sequencesDNAbin <- read.dna("/home/josephcrispell/Desktop/Research/RepublicOfIreland/Fastqs/vcfFiles/sequences_Prox-10_11-07-2018.fasta",
                            format = "fasta", skip=1) # skip first line - I added this line into FASTA: nSequences length

# Build the distance matrix
distanceMatrix <- dist.dna(sequencesDNAbin, model="JC69")

# Build a neighbour joining tree - an initial tree
njTree <- nj(distanceMatrix)

# Parse the isolate IDs
njTree$tip.label <- parseIDs(njTree$tip.label)

# Change the branch lengths to SNPs
njTree$edge.length <- njTree$edge.length * 632

# Remove NI isolates
njTreeWithoutNI <- drop.tip(njTree, c("Ref-1997", "182-MBovis", "161-MBovis"))

# Plot the phylogeny
pdf("/home/josephcrispell/Desktop/Research/RepublicOfIreland/Fastqs/vcfFiles/njTree_11-07-18.pdf")
par(mar=c(0,0,0,0))
plot.phylo(njTreeWithoutNI, label.offset=0.5, edge.color="grey", edge.width=2)

# Add a scale bar
axisLimits <- par("usr")
xPad <- 0.1* (axisLimits[2] - axisLimits[1])
yPad <- 0.1 * (axisLimits[4] - axisLimits[3])
lines(x=c(axisLimits[2] - xPad - 1, axisLimits[2] - xPad), y=c(axisLimits[3] + yPad, axisLimits[3] + yPad), lwd=2)
text(x=axisLimits[2] - xPad - 0.5, y=axisLimits[3] + (yPad*1.2), 
     labels="1 SNP")
dev.off()

#############
# FUNCTIONS #
#############

parseIDs <- function(tipLabels){
  
  newLabels <- c()
  for(i in 1:length(tipLabels)){
    
    newLabels[i] <- strsplit(tipLabels[i], split="_")[[1]][1]
  }
  
  return(newLabels)
}
