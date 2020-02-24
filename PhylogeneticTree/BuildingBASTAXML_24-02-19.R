#### Preparation ####

# Load the required libraries
library(ape)

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in FASTA and metadata ####

# Read in the FASTA file
fastaFile <- file.path(path, "example.fasta")
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)

# Create a metadata file


#### Select sequences ####

#### Build BASTA xml(s) ####

#### LOAD FUNCTIONS ####

