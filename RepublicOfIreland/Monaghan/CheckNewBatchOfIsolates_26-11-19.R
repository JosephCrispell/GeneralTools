# Set the path
path <- file.path("storage", "Research", "RepublicOfIreland", "Mbovis", "Monaghan")

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Read in the new batch of isolate IDs
file <- file.path(path, "IsolateIDs_26-11-19.txt")
isolates <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
isolates <- parseIsolateIDs(isolates$Isolate.ID)

# Check for duplicates
length(unique(isolates)) == length(isolates)

# Read in the Monaghan
file <- file.path(path, "SampleInformation_26-11-19.csv")
sampleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)
rownames(sampleInfo) <- sampleInfo$Aliquot

# Check if missing sample information for any isolate IDs
found <- sampleInfo[isolates[isolates %in% sampleInfo$Aliquot], ]
table(found$Species)
missing <- isolates[isolates %in% sampleInfo$Aliquot == FALSE]

# Read in the previous sample information table
file <- file.path(path, "SampleInformation_04-07-19.csv")
sampleInfoPrevious <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Check if any IDs in the previous sample information
table(isolates %in% sampleInfoPrevious$Aliquot)

#### FUNCTIONS ####

parseIsolateIDs <- function(ids){
  
  # Initialise a vector to stored the parsed IDs
  parsedIDs <- c()
  
  # Examine each of the isolate IDs
  for(id in ids){
    
    # Add a "TB" prefix
    id <- paste0("TB", id)
    
    # Split the ID into its parts
    parts <- strsplit(id, split="-")[[1]]
    
    # Pad second part with zeros to make length 6
    parts[2] <- paste0(paste(rep(0, 6-nchar(parts[2])), collapse=""), parts[2])
    
    # Create the complete new parsed label
    parsedIDs[length(parsedIDs) + 1] <- paste0(parts[1], "-", parts[2])
  }
  
  return(parsedIDs)
}