#### Preparation ####

# Load libraries
library(ape) # Build, plot and write NJ tree

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Create a path variable
path <- "/home/josephcrispell/Desktop/Research/Homoplasy/INDELS/"

#### Read in the data ####

# Read in the presence absence table
presenceAbsence <- read.table(paste0(path, "indel_sites.csv"), header=TRUE, sep=",", stringsAsFactors=FALSE)

# Get the transpose of the presence absence table
starts <- presenceAbsence$start
ends <- presenceAbsence$end
presenceAbsence <- t(presenceAbsence[, -c(1,2)])

#### Build a phylogeny ####

# Construct distance matrix
#   binary: distance is the proportion of bits in which only one is on amongst those in which at least one is on.
# dist calculate distance between rows
distances <- dist(presenceAbsence, method="binary")

# Build a neighbour joining tree
njTree <- nj(distances)

# Plot the phylogeny
plot.phylo(njTree)

# Write the phylogeny to file
write.tree(njTree, file=paste0(path, "neighbourJoining_", date, ".tree"))

#### FUNCTIONS ####