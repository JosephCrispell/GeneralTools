#### Load libraries ####

library(ape)

#### Read in the FASTA file ####

# Set the path variable
path <- "/home/josephcrispell/Desktop/"

# Read in the FASTA file
fastaFile <- paste0(path, "example_08-11-18.fasta")
sequences<- read.dna(fastaFile, format = "fasta")

# Note the number of sites in the FASTA file
nSitesInFasta <- length(as.list(sequences)[[1]])

#### Build a phylogeny ####

# Build the distance matrix
distanceMatrix <- dist.dna(sequences, model="JC69")

# Build a neighbour joining tree - an initial tree
njTree <- nj(distanceMatrix)

#### Plot the phylogeny ####

# Plot the tree
plot.phylo(njTree)

# Add node labels
nodelabels()

#### Calculate distance to MRCA ####

# Calculate a patristic distance matrix
patristicDistances <- dist.nodes(njTree)

# Pick an MRCA node - node with clade below
mrca <- 134 # Index for nodes in tips is 1:nTips and then (nTips+1):(nTips+1+tree$Nnode)

# Pick a tip
tip <- which(njTree$tip.label == "173")

# Note the number of tips
nTips <- length(njTree$tip.label)

# Calculate the distance in SNPs between the tip and its MRCA
patristicDistance <- patristicDistances[mrca, tip]

# Convert patristic distance to number of SNPs 
# NOTEs:
# - There won't be perfect agreement between phylogeny and sequences so SNP distances won't be whole numbers
# - As far as I am aware removing N sites is NOT appropriate. You could consider removing uninformative constant sites (same across all sequences)
nSNPs <- patristicDistance * nSitesInFasta
