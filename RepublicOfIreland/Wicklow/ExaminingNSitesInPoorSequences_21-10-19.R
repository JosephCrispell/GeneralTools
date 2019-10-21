# Load libraries
library(ape)
library(pathogenGenomicsWorkshop)

# Set the path 
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/Wicklow/vcfFiles/PreviousMergings/"

# Read in the FASTA file
fastaFile <- paste0(path, "sequences_Prox-10_19-03-2019.fasta")
sequences <- read.dna(fastaFile, format = "fasta", as.character=TRUE)

# Remove the sites that are Ns in 68-Mbovis and 44-2623
poor <- c("68-MBovis_31.vcf.gz", "44-2623_17.vcf.gz")
table(sequences["68-MBovis_31.vcf.gz", ]) / ncol(sequences)
table(sequences["44-2623_17.vcf.gz", ]) / ncol(sequences)

nSites <- unique(which(sequences["68-MBovis_31.vcf.gz", ] == "n" & sequences["44-2623_17.vcf.gz", ] == "n"))

sequences <- sequences[, -nSites]

distances <- dist.dna(as.DNAbin(sequences), model="JC69")

tree <- nj(distances)

tree <- drop.tip(tree, c("161-MBovis_2.vcf.gz", "182-MBovis_3.vcf.gz", "Ref-1997"))

plot(tree)
