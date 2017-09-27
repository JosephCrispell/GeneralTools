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

# Replace Clostroides difficile - NCBI changed name
bacteria[["Clostridioides difficile"]] <- NULL
bacteria[["Clostridium difficile"]] <- 0

bacteria <- countNumberOfTimesBacteriaAppearInTitles(bacteria,
                                                     popResults)

########################################
# Select the most represented bacteria #
########################################

# Get the count distribution
distribution <- getCountDistribution(bacteria)

# Note the top bacteria and their counts
topValues <- distribution[
  order(distribution, decreasing=TRUE)][1:10]
topNames <- names(bacteria)[
  order(distribution, decreasing=TRUE)][1:10]

##################
# Produce figure #
##################

# Produce an output figure
file <- paste(path, "CountsOfBacteriaInPoPResults_20-09-17.pdf", sep="")
pdf(file)

hist(distribution[distribution > 0], las=1, breaks=50,
     xlab="Number of hits", xlim=c(1,max(distribution)),
     main="Number of WGS research articles for bacteria")
legend("topright", legend=paste(topNames, " (", topValues, ")", sep=""),
       title="Top 10:", bty="n", cex=0.75)

dev.off()

#######################################
# Try and find bacteria in PoP titles #
#######################################

bacteria <- identifyPossibleScientificNames(popResults)

###############################################################
# Count how many times each bacteria appear in article titles #
###############################################################

bacteria <- countNumberOfTimesBacteriaAppearInTitles(bacteria,
                                                     popResults)

########################################
# Select the most represented bacteria #
########################################

# Two spellings of Klebsiella pneumoniae found
bacteria[["Klebsiella pneumoniae"]] <- bacteria[["Klebsiella pneumonia"]]
bacteria[["Klebsiella pneumonia"]] <- NULL

# Get the count distribution
distribution <- getCountDistribution(bacteria)

# Note the top bacteria and their counts
topValues <- distribution[
  order(distribution, decreasing=TRUE)][1:10]
topNames <- names(bacteria)[
  order(distribution, decreasing=TRUE)][1:10]

##################
# Produce figure #
##################


# Produce an output figure
file <- paste(path, "CountsOfBacteriaInPoPResults_FOUND_20-09-17.pdf", sep="")
pdf(file)

hist(distribution[distribution > 2], las=1, breaks=50,
     xlab="Number of hits", xlim=c(2,max(distribution)),
     main="Number of WGS research articles for bacteria")
legend("topright", legend=paste(topNames, " (", topValues, ")", sep=""),
       title="Top 10:", bty="n", cex=0.75)

dev.off()



#############
# FUNCTIONS #
#############

identifyPossibleScientificNames <- function(popResults){
  
  # Try and find bacteria names
  bacteria <- list()
  for(row in 1:nrow(popResults)){
    
    # Get title
    title <- popResults[row, "Title"]
    
    # Remove any punctuation from title
    title <- gsub("[[:punct:]]", " ", title)
    
    # Split it into words
    words <- strsplit(title, split=" ")[[1]]
    
    # Skip titles with less than two words
    if(length(words) == 1){
      next
    }
    
    # Try and find scientific name in title words
    for(i in 2:length(words)){
      
      # Get the first letters of words
      firstLetterOfFirstWord <- substr(words[i-1], start=1, stop=1)
      FirstLetterOfSecondWord <- substr(words[i], start=1, stop=1)
      
      # Check if current two words might be a scientific name
      if(nchar(words[i-1]) > 3 && nchar(words[i]) > 3 &&
         firstLetterOfFirstWord == toupper(firstLetterOfFirstWord) &&
         FirstLetterOfSecondWord == tolower(FirstLetterOfSecondWord)){
        
        # Build the name and store it
        name <- paste(words[i-1], words[i])
        if(is.null(bacteria[[name]]) == TRUE){
          bacteria[[name]] <- 0
        }
      }
    }
  }
  
  return(bacteria)
}

getCountDistribution <- function(bacteria){
  
  names <- names(bacteria)
  
  distribution <- c()
  for(i in 1:length(names)){
    distribution[i] <- bacteria[[names[i]]]
  }
  
  return(distribution)
}

countNumberOfTimesBacteriaAppearInTitles <- function(bacteria,
                                                     popResults){
  
  names <- names(bacteria)
  
  for(row in 1:nrow(popResults)){
    
    for(i in 1:length(names)){
      
      if(grepl(pattern=names[i], x=popResults[row, "Title"]) == TRUE){
        bacteria[[names[i]]] <- bacteria[[names[i]]] + 1
      }
    }
    
    if(row %% 100 == 0){
      cat(".")
    }
  }
  cat("\n")
  
  return(bacteria)
}

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