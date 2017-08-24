#################
# Load packages #
#################

library(seqinr)

#####################
# Set path variable #
#####################

# Set the path
path <- "/Users/josephcrisp1/Desktop/FilesToWorkOn/"

#####################
# Read in sequences #
#####################

# First sequence
file <- paste(path, "NC_002945.3_AF2122-97.fasta", sep="")
ref1Sequence <- as.vector(read.fasta(file)[[1]])

# Second sequence
file <- paste(path, "LT708304.1_AF2122-97.fasta", sep="")
ref2Sequence <- as.vector(read.fasta(file)[[1]])

################################
# Read in the annotation files #
################################

# Annotations for first sequence
file <- paste(path, "NC_002945_annotation.gff.txt", sep="")
annotations1 <- read.table(file, header=FALSE, skip=5, sep="\t", stringsAsFactors=FALSE)

# Annotations for the second sequence
file <- paste(path, "Garnier-MaloneTransfer.LT708304.1.Report.gff3", sep="")
annotations2 <- read.table(file, header=FALSE, skip=2, sep="\t", stringsAsFactors=FALSE)

########################################################################
# Find the locus tags - and their coordinates in each annotation table #
########################################################################

# Create a list to store the locus tags
locusCoordinates <- list()

# Examine annoation files
locusCoordinates <- findLocusInformation(annotations1, locusCoordinates, "Annotations1")
locusCoordinates <- findLocusInformation(annotations2, locusCoordinates, "Annotations2")

# Get info for locuses found in both annotation files
output <- findLociPresentInBoth(locusCoordinates)
lociFoundInBoth <- output[["Both"]]
lociFoundInSingle <- output[["Single"]]

###############################
# Compare the locus sequences #
###############################

# Note loci where either they are a different length or they differ genetically
problemTags <- compareLocusSequencesInReferences(ref1Sequence, ref2Sequence, locusCoordinates,
                                                 lociFoundInBoth)

#############
# FUNCTIONS #
#############

findLociOnlyFoundInOneAnnotation <- function()

compareLocusSequencesInReferences <- function(ref1Sequence, ref2Sequence, locusCoordinates,
                                              lociFoundInBoth){
  
  problemTags <- list()
  
  for(locusTag in lociFoundInBoth){
    
    # Get the sequence for the current locus from each sequence
    coordinates <- locusCoordinates[[locusTag]]
    sequence1 <- ref1Sequence[coordinates[["Annotations1"]][1]:coordinates[["Annotations1"]][2]]
    sequence2 <- ref2Sequence[coordinates[["Annotations2"]][1]:coordinates[["Annotations2"]][2]]
    
    # Compare the sequences
    geneticDist <- geneticDistance(sequence1, sequence2, locusTag)
    
    if(geneticDist != 0 && geneticDist != -1){
      cat(paste("Sequences weren't exactly the same!\t", locusTag, "\tdistance = ", geneticDist,
                "\n", sep=""))
      problemTags[[locusTag]] <- geneticDist
    }else if(geneticDist == -1){
      problemTags[[locusTag]] <- c(length(sequence1), length(sequence2), 
                                   abs(length(sequence1) - length(sequence2)))
    }
  }
  
  return(problemTags)
}

geneticDistance <- function(a, b, locusTag){
  
  dist <- 0
  
  if(length(a) != length(b)){
    cat(paste("Error! Loci aren't the same length\t", locusTag, "\n", sep=""))
    dist <- -1
  }else{
    
    for(i in 1:length(a)){
      
      if(a[i] != b[i]){
        dist <- dist + 1
      }
    }
  }
  
  return(dist)
}

findLociPresentInBoth <- function(locusCoordinates){
  
  # Get a list of all the locus tags
  locusTags <- names(locusCoordinates)
  
  # Initialise a list to store locus coordinates for loci found in both annotations
  lociFoundInBoth <- c()
  lociFoundInOne <- list()
  
  # Find which ones weren't present in both
  for(locusTag in locusTags){
    
    if(length(names(locusCoordinates[[locusTag]])) == 2){
      
      lociFoundInBoth[length(lociFoundInBoth) + 1] <- locusTag
    }else{
      lociFoundInOne[[locusTag]] <- names(locusCoordinates[[locusTag]])[1]
    }
  }
  
  output <- list(
    "Both"=lociFoundInBoth,
    "Single"=lociFoundInOne
  )
  
  return(output)
}

findLocusInformation <- function(annotationTable, locusCoordinates, tableTag){
  
  for(row in 1:nrow(annotationTable)){
    
    # Get the locus information from the current row
    locusInfo <- annotationTable[row, ncol(annotationTable)]
    
    # Split the locus information
    parts <- strsplit(locusInfo, split=";")[[1]]
    
    # Examine the information available for the locus and look for the locus tag
    locusTag <- NA
    for(i in 1:length(parts)){
      
      # Look for the locus tag
      if(grepl(x=parts[i], pattern="locus_tag") == TRUE){
        
        # Get the locus tag id
        locusTag <- strsplit(parts[i], split="=")[[1]][2]
        
        # Store the coordinates for the locus
        if(is.null(locusCoordinates[[locusTag]]) == TRUE){
          
          locusCoordinates[[locusTag]][[tableTag]] <- 
            c(annotationTable[row, 4], annotationTable[row, 5])
          
        }else if(tableTag %in% names(locusCoordinates[[locusTag]]) == FALSE){
          
          locusCoordinates[[locusTag]][[tableTag]] <- 
            c(annotationTable[row, 4], annotationTable[row, 5])
        }
        break
      }
    }
  }
  
  return(locusCoordinates)
}

