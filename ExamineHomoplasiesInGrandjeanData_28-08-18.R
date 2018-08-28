#### Load Libraries ####
library(devtools)
library(homoplasyFinder)

#### Note FASTA and NEWICK file locations ####

# Create path variable
path <- "/home/josephcrispell/Desktop/Research/Granjean2017ConvergentEvolution/"

# Note FASTA file
fastaFile <- paste0(path, "FASTQs/vcfFiles/sequences_Prox-0_27-08-2018.fasta")

# Note NEWICK file - NOTE: removed ">" from tip labels
newickFile <- paste0(path, "RAxML_27-08-18/RAxML_bipartitions.RaxML-R_27-08-18_rmGT")

#### Run HomoplasyFinder ####

# Run homoplasyFinder
inconsistentSites <- runHomoplasyFinderInJava(treeFile=newickFile, fastaFile=fastaFile, path=path)

# Read in consistency index report
results <- read.table(paste0(path, "HomoplasyFinderOutput/consistencyIndexReport_28-08-18.txt"), header=TRUE, sep="\t")

#### Identify genome positions of inconsistent sites ####

# Read in the FASTA positions file
fastaPositions <- read.table(paste0(path, "FASTQs/vcfFiles/fastaPositions_Prox-0_27-08-2018.txt"), 
                             header=TRUE)

# Add genome position to results table
results$GenomePosition <- fastaPositions[results$Position, 1]

#### Read in table of homoplasies found by Grandjean et al. ####

# Read in table of most homoplastic polymorphisms
mostHomoplastic <- read.table(paste0(path, "SupplementaryInformation/MostHomoplasticPolymorphisms_28-08-18.csv"), 
                              header=TRUE, sep=",")

