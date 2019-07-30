# Set the path
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/"

# Read in the full Wicklow isolate table
file <- paste0(path, "IsolateSpeciesAndYear_26-04-19.csv")
wicklow <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
wicklow$ID <- getIDFromWicklowLabels(wicklow$Aliquot)

# Read in the latest batch
file <- paste0(path, "IsolateBatch_11-06-19.csv")
batch <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
batch$ID <- editLabels(batch$Isolate.ID)

# Compare IDs
idsInWicklow <- batch$ID[which(batch$ID %in% wicklow$ID)]
length(idsInWicklow)

subset <- wicklow[wicklow$ID %in% batch$ID, ]

tipInfo$ID2 <- getIDFromWicklowLabels(tipInfo$Aliquot)

#### FUNCTIONS ####

getIDFromWicklowLabels <- function(labels){
  
  # Create a vectore to store the output
  output <- c()
  
  # Exmaine each of the IDs
  for(index in seq_along(labels)){
    
    # Get the second part of the label and store
    output[index] <- strsplit(labels[index], split="-")[[1]][2]
  }
  
  return(output)
}

editLabels <- function(ids){
  
  # Initialise a vector to store the edited IDs
  output <- c()
  
  # Examine each ID
  for(index in seq_along(ids)){
    
    # Pad with zeros to make it 8 characters
    output[index] <- paste0(paste(rep(0, 6 - nchar(ids[index])), collapse=""), ids[index])
  }
  
  return(output)
}