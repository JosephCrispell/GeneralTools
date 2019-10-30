#### Preparation ####

# Load the libraries
library(ape)
library(phytools) # to.matrix

# Set the path
path <- file.path("~", "Desktop")

#### Read in the Average Nucleotide Identity ####

# Read in the Average Nucleotide Identify table
aniFile <- file.path(path, "ani_matrix.csv")
averageNucleotideIdentity <- read.table(aniFile, header=TRUE, sep=",", stringsAsFactors=FALSE, row.names=1)
species <- averageNucleotideIdentity["species"]
averageNucleotideIdentity <- averageNucleotideIdentity[, -which(colnames(averageNucleotideIdentity) == "species")]


#### Read in the region hits data ####

# Read in the region hits
regionsFile <- file.path(path, "rd900_region_hits.csv")
regions <- read.table(regionsFile, header=TRUE, sep=",", stringsAsFactors=TRUE, row.names=4)


#### Construct the phylogeny ####

# Create the distance matrix
distances <- as.matrix(1 - (averageNucleotideIdentity/100))

# Build a neighbour joining tree
njTree <- nj(distances)

# Root the tree
rootedTree <- root(njTree,
                   outgroup="GCA_001287165", # Set the isolate to root the tree on - can also use node
                   resolve.root=TRUE) # Ensures new root will be bifurcating

# Remove any non-dichotomous branching nodes (nodes with more than two daughters)
rootedTree <- multi2di(rootedTree, random=TRUE)

# Make branches with lengths of 0 very very small
rootedTree$edge.length[rootedTree$edge.length <= 0] <- 0.000001

#### Get the tip information ####

# Create a tipInfo table
tipInfo <- regions[, c("pro_pknh1", "pro_pknh2", "pknh2_sensor")]
tipInfo <- tipInfo[rootedTree$tip.label, ]

# Add states to the tip info table
tipInfo[tipInfo > 0] <- 1
tipInfo$states <- as.factor(paste(tipInfo$pro_pknh1, tipInfo$pro_pknh2, tipInfo$pknh2_sensor, sep="-"))

# Add the species
tipInfo$species <- species[rootedTree$tip.label, "species"]

#### Fit the ancestral character states ####

# Estimate ancestral states under a ER model
ancestralStateFitting <- ace(tipInfo$states, rootedTree, 
                             model="ARD", # "ER" Equal rates; "ARD" All rates different
                             type="discrete") # Discrete states

#### Plot the phylogeny and estimated states ####

# Set the plotting margins: c(bottom, left, top, right)
par(mar=c(0, 0, 0, 1))

# Assign colours to the states
stateColours <- setNames(c('red','blue','green2','orange','purple','yellow'), levels(tipInfo$states))

# Plot the phylogeny
plot.phylo(rootedTree, type="phylogram", # Set the shape of the tree to be a phylogram
           show.tip.label=FALSE) # Don't plot the tip labels

# Add pie charts at each internal node
nodelabels(node=1:rootedTree$Nnode+Ntip(rootedTree), # Provide an array with all the internal node IDs (nodes numbered from nTips onwards...)
           pie=round(ancestralStateFitting$lik.anc, 3), # Assign pie charts to each node - based on the outputs from ancestral state estimation
           piecol=stateColours, # Set the colours to be used in the pie charts
           cex=0.3) # Set the size of the pie charts

# Add circles to indicate the states of the tips
tiplabels(pie=to.matrix(tipInfo$states, levels(tipInfo$states)), # Provide the state of each tip
          piecol=stateColours, # Provide the colour assigned to each state
          cex=0.2) # Set the size of the tip circle

# Add a legend
axisLimits <- par()$usr
add.simmap.legend(colors=stateColours, prompt=FALSE,
                  x=0.75*axisLimits[2],
                  y=0.9*axisLimits[4], fsize=1.2)

# Add species labels
tiplabels(text=tipInfo$species, offset=0.00005, adj=0, frame='n', cex=0.9, xpd=TRUE)

