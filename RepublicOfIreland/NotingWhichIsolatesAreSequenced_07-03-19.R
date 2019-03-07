#### Read in the sample information ####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/RepublicOfIreland/Mbovis/"

# Load in the table with the CVRL Aliquot codes
file <- paste0(path, "all_lab_samples.csv")
allCodes <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Create a parse set of aliquot codes to match those in the to be sequenced table
parsedCodes <- parseCodes(allCodes$Aliquot)

# Read in the table detailing the sequenced isolates
file <- paste0(path, "Mbovis_CattleSamplingInfo_17-07-18.tsv")
sequenced <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read in the table detailing the isolates to be sequenced
file <- paste0(path, "Feb2019SequencingFinal.csv")
toBeSequenced <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#### Note which isolates have been/are being sequenced ####

# Add a column to the table with all the isolate codes
allCodes$Progress <- "Require heat-killed cells"

# Note the sequenced isolates
allCodes[allCodes$Aliquot %in% sequenced$Aliquot, "Progress"] <- "Sequenced"

# Note the isolates just about to be sequenced - CAN'T GET FORMAT OF IDS TO MATCH
allCodes[parsedCodes %in% toBeSequenced$Sample.isolate.ID, "Progress"] <- "Heat-killed cells received and DNA prepared for sequencing"

#### FUNCTIONS ####

parseCodes <- function(codes){
  
  # Initialise a vector to store the output
  output <- c()
  
  # Examine each code
  for(i in seq_along(codes)){
    
    # Remove the "TB"
    output[i] <- substr(codes[i], 3, nchar(codes[i]))
    
    # Get rid of "-00"
    output[i] <- gsub(output[i], pattern="-00", replacement="-")
    
    # Get rid of "-0"
    output[i] <- gsub(output[i], pattern="-0", replacement="-")
  }
  
  return(output)
}