##################
# Get input data #
##################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/LiteratureReview/"

# Read in the Publish Or Perish table
file <- paste(path, "PoP_Results_2007-2017_08-09-17.csv", sep="")
popResults <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the NCBI list of bacteria with complete genomes
file <- paste(path, "NCBI_Bacteria_CompleteGenomes.txt", sep="")
ncbiBacteria <- readTable(file, sep="\t")

###########################################################################
# Get a list of bacteria names to search for in Publish or Perish results #
###########################################################################

bacteria <- getListOfBacteria(ncbiBacteria)

###############################################################
# Count how many times each bacteria appear in article titles #
###############################################################

names <- names(bacteria)

for(row in 1:nrow(popResults)){
  
  for(i in 1:length(names)){
    
    #parts <- strsplit(names[i], split = " ")[[1]]
    
    if(grepl(pattern=names[i], x=popResults[row, "Title"]) == TRUE){
      bacteria[[names[i]]] <- bacteria[[names[i]]] + 1
    }
    
    #if(grepl(pattern=parts[2], x=popResults[row, "Title"]) == TRUE){
    #  bacteria[[names[i]]] <- bacteria[[names[i]]] + 1
    #}
  }
  
  if(row %% 100 == 0){
    cat(".")
  }
}
cat("\n")

########################################
# Select the most represented bacteria #
########################################

distribution <- c()
for(i in 1:length(names)){
  distribution[i] <- bacteria[[names[i]]]
}

hist(distribution[distribution > 0], las=1, breaks=50,
     xlab="Number of hits", main="")
head(names[order(distribution, decreasing=TRUE)], n=20)


#############
# FUNCTIONS #
#############

getListOfBacteria <- function(ncbiBacteria){
  
  # Initialise a list to store the bacteria names
  bacteria <- list()
  
  # Examine each of the bacteria listed in the NCBI table
  for(row in 1:nrow(ncbiBacteria)){
    
    # Get the first two parts of the bacteria's name
    parts <- strsplit(ncbiBacteria[row, "#Organism/Name"], split = " ")[[1]]
    name <- paste(parts[1], parts[2], sep=" ")
    
    # Remove any punctuation in the names
    name <- gsub("[[:punct:]]", "", name)
    
    # Check if already exists in list - if not then add it
    if(is.null(bacteria[[name]]) == TRUE){
      bacteria[[name]] <- 0
    }
  }
  
  return(bacteria)
}

readTable <- function(file, sep){
  
  # Open the file and store all the lines
  connection <- file(file, open="r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Split the array of file lines into an array
  table <- data.frame(do.call(rbind, strsplit(fileLines[2:length(fileLines)], split=sep)))
  
  # Set the column names
  colnames(table) <- strsplit(fileLines[1], split=sep)[[1]]
  
  # Convert each column to character vector
  for(col in 1:ncol(table)){
    table[, col] <- as.character(table[, col])
  }
  
  return(table)
}