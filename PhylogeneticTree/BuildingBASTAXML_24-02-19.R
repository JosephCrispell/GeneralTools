#### Preparation ####

# Load the required libraries
library(ape) # Reading in FASTA
library(lubridate) # Converting to decimal dates

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in FASTA, metadata and constant site counts ####

# Read in the FASTA file
fastaFile <- file.path(path, "example.fasta")
sequences <- read.dna(fastaFile, format="fasta", as.character=TRUE)

# Create a metadata file
generateRandomMetadata(nrow(sequences), path)

# Read in the metadata
metadataFile <- file.path(path, "metadata.csv")
metadata <- read.table(metadataFile, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Convert sampling dates to decimal dates
metadata$DecimalDate <- decimal_date(as.Date(metadata$Date))

# Read in the constant site counts
# - Could estimate these: (genomeSize - nSites) split between A, C, G, T with 65% bias to G and C

#### Build BASTA xml(s) ####

# Convert the sampling dates

#### LOAD FUNCTIONS ####

generateConstantSiteCounts <- function(nSites, genomeSize=4345492, nucleotideProbs=c(0.25, 0.25, 0.25, 0.25))

generateRandomMetadata <- function(nSequences, path){
  
  # Generate random data for sequences
  metadata <- data.frame("ID"=rownames(sequences), 
                         "Date"=sample(seq(from=as.Date("01-01-2005", format="%d-%m-%Y"),
                                           to=as.Date("01-12-2015", format="%d-%m-%Y"), by=1),
                                       size=nSequences, replace=TRUE),
                         "Species"=sample(c("badger", "cow"), size=nSequences, replace=TRUE))
  
  # Write the data to file
  metadataFile <- file.path(path, "metadata.csv")
  write.table(metadata, file=metadataFile, quote=FALSE, sep=",", row.names=FALSE)
}
