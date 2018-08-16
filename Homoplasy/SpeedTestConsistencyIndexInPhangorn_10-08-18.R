detach("package:phangorn", unload=TRUE)
detach("package:geiger", unload=TRUE)
detach("package:ape", unload=TRUE)

start.time <- Sys.time()

library(ape)
library(phangorn)

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Homoplasy/"

# Read in the tree file
tree <- read.tree(paste(path, "example-AFTER_10-08-18.tree", sep=""))

# Read in the sequence file
sequences <- read.phyDat(paste(path, "example_10-08-18.fasta", sep=""),format="fas")

# Calculate the consistency index
ci <- CI(tree, sequences, sitewise=TRUE)

# Identify inconsistent sites
inconsistent <- (which(!is.nan(ci) & ci < 1))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

detach("package:phangorn", unload=TRUE)
detach("package:geiger", unload=TRUE)
detach("package:ape", unload=TRUE)
inconsistent
