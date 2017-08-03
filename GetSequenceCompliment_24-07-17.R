# Set path
path <- "C:/Users/Joseph Crisp/Desktop/"

# Read in sequence
#sequence <- readSequence(paste(path, "RD1_sequence_AF2122-97.txt", sep=""))

sequence <- "GTCGTCAGACCCAAAA"

# Convert sequence to characters
nucleotides <- strsplit(sequence, "")[[1]]

# Reverse the nucleotides
reverse <- rev(nucleotides)

# Get the nucleotide sequence compliment
reverseCompliment <- getCompliment(reverse)

# Print resulting sequence
print(paste(nucleotides, collapse=""))
print(paste(reverseCompliment, collapse=""))

############
# FUNCTION #
############

getCompliment <- function(nucleotides){
  
  # Create a list for the nucleotide compliments
  complimentNucleotide <- list(
    'A'='T',
    'C'='G',
    'G'='C',
    'T'='A'
  )
  
  # Initialise an array to store each nucleotide's compliment
  compliment <- c()
  
  # Get the compliment for each nucleotide
  for(i in 1:length(nucleotides)){
    compliment[i] <- complimentNucleotide[[nucleotides[i]]]
  }
  
  return(compliment)
}

readSequence <- function(fileName){
  
  # Read lines from file
  connection <- file(fileName, open="r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Return first line - sequence
  return(fileLines[1])
}