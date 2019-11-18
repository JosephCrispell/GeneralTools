#### Preparation ####

# Load libraries
library(ape)

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the path
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/Wicklow/"

#### Read in the data ####

# Read in table that links original sequence ID to aliquot IDs
file <- paste0(path, "Mbovis_SamplingInfo_17-07-18.tsv")
linkTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Read in the isolate metadata
file <- paste0(path, "IsolateSpeciesAndYear_26-04-19.csv")
metadata <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")

#### Read in the phylogeny ####

# Read in the phylogeny
file <- paste0(path, "RAxML_21-10-19/", "RAxML_bipartitions.RaxML-R_21-10-19")
tree <- read.tree(file)

# Remove NI isolates and Reference
tree <- drop.tip(tree, tree$tip.label[grepl(tree$tip.label, pattern=">Ref-1997|>182-MBovis|>161-MBovis")])

# Parse the tip labels
ids <- editTipLabels(tree$tip.label)

#### Get isolate information ####

# Get the sampling date and species for each isolate ID
isolateInfo <- getIsolateInfo(ids, metadata, linkTable)

#### Create the biosample ####

# Read in the template file
templateFile <- paste0(path, "SRA_submission_18-11-19/", "Biosample_template_18-11-19.tsv")
template <- read.table(templateFile, header=TRUE, check.names=FALSE, skip=11, sep="\t", stringsAsFactors=FALSE)

# Create a table to hold the badger and cattle info
output <- as.data.frame(matrix(NA, nrow=nrow(isolateInfo), ncol=ncol(template)), stringsAsFactors=FALSE)
colnames(output) <- colnames(template)

# Fill the table with the cattle and badger sequence data
output$`*sample_name` <- isolateInfo$ID
output$isolate <- isolateInfo$ID
output$`*organism` <- rep("Mycobacterium tuberculosis variant bovis", nrow(isolateInfo))
output$`*collection_date` <- isolateInfo$Date
output$`*geo_loc_name` <- rep("Ireland: County Wicklow", nrow(isolateInfo))
output$`*sample_type` <- rep("Culture", nrow(isolateInfo))
output$host <- toupper(isolateInfo$Species)

# Write the output table to file
outputFile <- paste0(path, "SRA_submission_18-11-19/", "Biosample_Wicklow_18-11-19.tsv")
write.table(output, file=outputFile, quote=FALSE, sep="\t", row.names=FALSE, na="")

#### FUNCTIONS ####

getIsolateInfo <- function(isolates, metadata, linkTable){
  
  # Initialise a dataframe to store the tip information
  isolateInfo <- data.frame(ID=isolates, Species=NA, Date=NA, stringsAsFactors=FALSE)
  
  # Examine each of the tips
  for(index in seq_along(isolates)){
    
    # Initialise variables to store the tip's information
    aliquotCode <- NA
    species <- NA
    date <- NA

    # Check if the current tip is associated with the original dataset
    if(grepl(isolates[index], pattern="-MBovis")){
      
      # Get the sequence number from the curren tip label
      sequenceNumber <- strsplit(isolates[index], split="-")[[1]][1]
      
      # Find the row in the link table
      row <- which(linkTable$Isolate.Code == sequenceNumber)
      
      # Get the current tips aliquot code
      if(length(row) != 0){
        aliquotCode <- linkTable[row, "Aliquot"]
      }else{
        cat(paste("Error for old batch. Couldn't find sequence number: ", sequenceNumber, " ", isolates[index], "\n"))
      }
      
      # Get information from most recent sequencing run
    }else{
      
      # Get the second part of the tip label - looks like an aliuot label without 00s
      aliquotCodePart <- strsplit(isolates[index], split="-")[[1]][2]
      
      # Find row that matches above part
      row <- which(grepl(metadata$Aliquot, pattern=aliquotCodePart))
      
      # Get the full aliquot code
      if(length(row) == 1){
        aliquotCode <- metadata[row, "Aliquot"]
      }else{
        cat(paste("Error for new batch. Couldn't find aliquot part: ", aliquotCodePart, " (found ", length(row), " matches)\n\n"))
      }
    }
    
    # If aliquot code available then get the tip species and sampling year
    if(is.na(aliquotCode) == FALSE){
      
      # Get the row in the metadata table for the current aliquot code
      row <- which(metadata$Aliquot == aliquotCode)
      
      # Note the species and convert multiple "Cow labels to single "Cow label
      species <- metadata[row, "Species"]
      if(species %in% c("Heifer", "Steer", "Calf", "Bull", "Bovine")){
        species <- "Cow"
      }
      
      # Note the sampling date
      if(metadata[row, "Received.Test.date"] != ""){
        date <- as.character(as.Date(metadata[row, "Received.Test.date"], format="%d/%m/%Y"))
      }
    }
    
    # Check species isn't nothing
    if(is.na(species) == FALSE && species == ""){
      species <- NA
    }

    # Store the current tips information
    isolateInfo[index, "Species"] <- species
    isolateInfo[index, "Date"] <- date
  }
  
  return(isolateInfo)
}

editTipLabels <- function(tipLabels){
  
  # Initialise a vector to store the new labels
  output <- c()
  
  # Examine each tip label
  for(label in tipLabels){
    
    # Check if starts with ">"
    if(grepl(label, pattern="^>")){
      
      # Remove the ">" prefix
      label <- substr(label, 2, nchar(label))
    }
    
    # Split the label and retain first part
    label <- strsplit(label, split="_")[[1]][1]
    
    # Store the new label
    output[length(output) + 1] <- label
  }
  
  return(output)
}
