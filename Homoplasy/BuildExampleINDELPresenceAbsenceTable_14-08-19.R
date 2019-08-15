
# Read in the FASTA file
fastaFile <- system.file("extdata", "example.fasta", package = "homoplasyFinder")
sequences <- readFASTA(fastaFile)

# Define some example regions
starts <-        c(35, 68, 92, 100, 150, 234, 270, 343, 399, 468)
ends <- starts + c(10, 15, 2,  45,  1,   25,  64,  12,  1,   500)

# Build the INDEL presence/absence table
presenceAbsence <- data.frame(matrix(sample(c(0,1), size=length(starts) * (length(sequences)+2), replace=TRUE),
                                     nrow=length(starts), ncol=length(sequences)+2), stringsAsFactors=FALSE)
colnames(presenceAbsence) <- c("start", "end", names(sequences))
presenceAbsence$start <- starts
presenceAbsence$end <- ends

# Take a quick look
presenceAbsence[, 1:5]

# Write to file
outputFile <- paste0(system.file("extdata", package = "homoplasyFinder"), "/presenceAbsence_INDELs.csv")
write.table(presenceAbsence, file=outputFile, sep=",", quote=FALSE, row.names=FALSE)

#### FUNCTIONS ####

readFASTA <- function(fastaFile, skip=1){
  
  # Open a connection to a file to read (open="r")
  connection <- file(fastaFile, open="r")
  
  # Get all lines from file and store in vector
  fileLines <- readLines(connection)
  
  # Close file connection
  close(connection)
  
  # Skip the first X lines
  fileLines <- fileLines[-skip]
  
  # Initialise a list to store the sequences
  sequences <- list()
  
  # Loop through each of the lines in file
  for(line in fileLines){
    
    # Check if found new sequence
    if(grepl(line, pattern="^>")){
      
      # Get the sequence name
      name <- substr(line, 2, nchar(line))
      
      # Start storing the sequence
      sequences[[name]] <- c()
      
      # Else if haven't found new sequence - continue storing previous sequence
    }else{
      sequences[[name]] <- c(sequences[[name]], strsplit(line, split="")[[1]])
    }
  }
  
  return(sequences)
}
