# Load the ape library
library(ape)

# Note the FASTA file name and path
fastaFile <- "/home/josephcrispell/Desktop/Research/Granjean2017ConvergentEvolution/FASTQs/vcfFiles/sequences_Prox-0_27-08-2018.fasta"

# Set the working directory
setwd("/home/josephcrispell/Desktop/Research/Granjean2017ConvergentEvolution/")

# Build analysis name
analysisName <- "RaxML-R_27-08-18"

# Set the input parameters for RAXML
model <- "GTRCAT" # No rate heterogenity
seeds <- sample(1:100000000, size=2, replace=FALSE) # For parsimony tree and boostrapping
nBootstraps <- 100
nThreads <- 10

# Build the command for RAXML
command <- paste("raxmlHPC-PTHREADS", # Note on WINDOWS replace with: /path/to/raxmlHPC-PTHREADS.exe
                 " -f a", # Algorithm: Rapid boostrap inference
                 " -N ", nBootstraps,
                 " -T ", nThreads,
                 " -m ", model, " -V", # -V means no rate heterogenity
                 " -p ", seeds[1], " -x ", seeds[2], # Parsimony and boostrapping seeds
                 " -n ", analysisName,
                 " -s ", fastaFile, sep="")

# Run RAXML
system(command, intern=TRUE)

# Get files in current working directory
files <- list.files()

# Select the tree file with BS support values
treeBSFile <- files[grepl(files, pattern=paste("RAxML_bipartitions[.]", analysisName, sep="")) == TRUE]

# Open the file
treeBS <- read.tree(treeBSFile)

treeBSFile
