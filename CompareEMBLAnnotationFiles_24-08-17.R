##################
# Load libraries #
##################

library(Biostrings)

##########################
# Read in the EMBL files #
##########################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Reference/"

# Original
# file <- paste(path, "TransferAnnotations/embl/", "AF2122-97_Garnier2003.embl", sep="")
# original <- readEMBLFile(file)
file <- paste(path, "UpdatedReference_Malone2017/",
              "LT708304.1_AF2122-97_Malone2017_rmPrefixToLocusTag.embl", sep="")
original <- readEMBLFile(file)

# Transferred
file <- paste(path, "TransferAnnotations/",
              "Garnier-MaloneTransfer_24-08-17.LT708304.1.final.embl", sep="")
transferred <- readEMBLFile(file)

########################
# Compare the features #
########################

############################################################################################
# Alignment algorithm options:                                                             #
#   global - Needlieman-Wunsch    local - Smith-Waterman    overlap - ends-free
#
#   global = align whole strings with end gap penalties                                    #
#   local = align string fragments                                                         #
#   overlap = align whole strings without end gap penalties                                #
#   global-local = align whole strings in pattern with consecutive subsequence of subject  #
#   local-global = align consecutive subsequence of pattern with whole strings in subject. #
############################################################################################


# Get common features
commonFeatures <- getCommonFeatures(original, transferred)

# Note the sequence of each feature
featureSequences <- getFeatureSequencesFromBothAnnotations(original, transferred)

# Compare the feature sequences between the annotation reports
substitutionMatrix <- nucleotideSubstitutionMatrix(match=1, mismatch=-1, 
                                                   baseOnly=TRUE)
alignmentInfo <- calculateAlignmentScores(commonFeatures, featureSequences,
                                          gapOpeningPenalty=0,
                                          gapExtensionPenalty=1,
                                          substitutionMatrix=substitutionMatrix,
                                          alignmentType="overlap")
alignmentScores <- alignmentInfo[["Scores"]]

#########################################
# Examine features with poor alignments #
#########################################

# Get a list of features without perfect alignment
threshold <- 0.95
poorFeatureAlignments <- alignmentScores[alignmentScores$Score < threshold, "Feature"]
length(poorFeatureAlignments)

# Examine each of the poorly aligned features
examineFeaturesWithPoorAlignments(poorFeatureAlignments, alignmentInfo)

#############
# FUNCTIONS #
#############

plotFeatureTypes <- function(features, featureInfo, main){
  
  # Get the types of the current features
  types <- c()
  for(i in 1:length(features)){
    types[i] <- featureInfo[[features[i]]]$type
  }
  
  # Make a bar plot
  barplot(table(types), las=2, cex.names=0.6, main=main)
}

examineFeaturesWithPoorAlignments <- function(poorFeatureAlignments, alignmentInfo){
  
  alignmentScores <- alignmentInfo[["Scores"]]
  
  for(feature in poorFeatureAlignments){
    cat("################################################################################################################\n")
    cat("################################################################################################################\n\n")
    alignment <- alignmentInfo[["Alignments"]][[feature]]
    print(alignment)
    cat(paste(score(alignment) / nchar(alignment), "\n\n\n\n"))
  }
  
  # Plot the scores
  hist(alignmentScores[alignmentScores$Score < threshold, "Score"], las=1, breaks=50,
       main="Features with poor alignments", xlab="Alignment Score")
}

calculateAlignmentScores <- function(commonFeatures, featureSequences,
                                     gapOpeningPenalty, gapExtensionPenalty,
                                     substitutionMatrix, alignmentType){
  
  # Initialise an output table to store the alignment scores
  alignmentScores <- data.frame(Feature=commonFeatures, stringsAsFactors=FALSE,
                                Score=rep(NA, length(commonFeatures)))
  
  alignments <- list()
  
  # Compare the sequences avalable for each feature
  for(i in 1:length(commonFeatures)){
    
    feature <- commonFeatures[i]
    
    alignment <- NULL
    
    # Compare the sequences for the current feature for the original and transferred
    if(nchar(featureSequences[[feature]]$Original) >= 
       nchar(featureSequences[[feature]]$Transferred)){
      alignment <- pairwiseAlignment(
        pattern=featureSequences[[feature]]$Original,
        subject=featureSequences[[feature]]$Transferred,
        scoreOnly=FALSE, type=alignmentType, 
        substitutionMatrix=substitutionMatrix,
        gapOpening=gapOpeningPenalty, gapExtension=gapExtensionPenalty)
    }else{
      alignment <- pairwiseAlignment(
        pattern=featureSequences[[feature]]$Transferred,
        subject=featureSequences[[feature]]$Original,
        scoreOnly=FALSE, type=alignmentType, 
        substitutionMatrix=substitutionMatrix,
        gapOpening=gapOpeningPenalty, gapExtension=gapExtensionPenalty)
    }
    
    # Store the alignment
    alignments[[feature]] <- alignment
    
    # Store the scores
    alignmentScores[i, "Score"] <- score(alignment) / nchar(alignment)
    
    # Progress
    if(i %% 1000 == 0){
      cat(paste("Finished alignment for", i, "of", length(commonFeatures), 
                "features.\n"))
    }
  }
  
  output <- list(
    "Scores"=alignmentScores,
    "Alignments"=alignments
  )
  
  return(output)
}

getFeatureSequencesFromBothAnnotations <- function(original, transferred){
  
  featureSequences <- list()
  for(feature in commonFeatures){
    
    # Store the sequence of the current feature in the original
    # and transferred annotations
    featureSequences[[feature]] <- list(
      "Original"=getFeatureSequence(feature, original),
      "Transferred"=getFeatureSequence(feature, transferred)
    )
  }
  
  return(featureSequences)
}

getFeatureSequence <- function(feature, annotationInfo){

  coordinates <- annotationInfo[["Features"]][[feature]]$coordinates
  
  # Get the corresponding sequence
  sequence <- c()
  for(i in seq(1, length(coordinates), 2)){
    
    # Check if complement or forward direction
    if(annotationInfo[["Features"]][[feature]]$direction == "FORWARD"){
      sequence <- c(sequence, annotationInfo[["Sequence"]][
        coordinates[i]:coordinates[i+1]])
      
    }else if(annotationInfo[["Features"]][[feature]]$direction == "COMPLEMENT"){
      sequence <- c(sequence, annotationInfo[["Sequence"]][
        coordinates[i+1]:coordinates[i]])
      
    }else{
      cat(paste("ERROR: Feature direction not recognised:", 
                annotationInfo[["Features"]][[feature]]$direction, "\n"))
    }
  }
  
  # Convert the sequence to a string
  sequence <- paste(sequence, collapse="")
  
  return(sequence)
}

getCommonFeatures <- function(original, transferred){
  
  # Get the common features
  commonFeatures <- intersect(names(original$Features),
                              names(transferred$Features))
  
  # Check features unique to either
  featuresUniqueToOriginal <- names(original$Features)[
    !(names(original$Features) %in% names(transferred$Features))]
  featuresUniqueToTransferred <- names(transferred$Features)[
    !(names(transferred$Features) %in% names(original$Features))]
  
  cat(paste("Found", length(commonFeatures), "common features.\n"))
  cat(paste("Found", length(featuresUniqueToOriginal), 
            "features unique to original.\n"))
  cat(paste("Found", length(featuresUniqueToTransferred), 
            "features unique to transferred.\n"))
  
  # Examine features unique to each annotation file
  if(length(featuresUniqueToOriginal) != 0){
    plotFeatureTypes(featuresUniqueToOriginal, original$Features, 
                     main="Types of Features Unique to Original")
  }
  if(length(featuresUniqueToTransferred) != 0){
    plotFeatureTypes(featuresUniqueToTransferred, transferred$Features, 
                     main="Types of Features Unique to Transferred")
  }
  
  return(commonFeatures)
}

readEMBLFile <- function(fileName){
  
  # Read in all the file lines
  connection <- file(fileName, open="r")
  fileLines <- readLines(connection)
  close(connection)
  
  # Initialise variables for recording feature information
  features <- list() # locus_tag, type, coordinates, direction
  foundFeature <- FALSE
  foundSequence <- FALSE
  sequenceParts <- c()
  locusTag <- NA
  geneTag <- NA
  
  # Parse the file lines
  for(i in 1:length(fileLines)){

    # Split the current line into its columns
    cols <- strsplit(fileLines[i], split="  +")[[1]]

    ### Parsing features ###
    
    # Note if found feature
    if(length(cols) == 3 && 
       grepl(cols[3], pattern="..") == TRUE &&
       cols[2] %in% c("gene", "CDS", "tRNA", "repeat_region",
                      "mobile_element", "rRNA", "misc_RNA")){
      
      # Store the information for the previous feature
      if(is.na(locusTag) == FALSE && grepl(locusTag, pattern="XXXX") == FALSE){
        features[[locusTag]] <- list(
          "locus_tag"=locusTag,
          "gene_tag"=geneTag,
          "type"=type,
          "direction"=direction,
          "coordinates"=coordinates)
      }else if(is.na(geneTag) == FALSE){

        cat(paste("Locus tag:", locusTag, "\tGene tag:", geneTag,
                  "\tType:", type, "\n"))
        
        features[[geneTag]] <- list(
          "locus_tag"=locusTag,
          "gene_tag"=geneTag,
          "type"=type,
          "direction"=direction,
          "coordinates"=coordinates)
      }
      
      # Reset the tag variables
      locusTag <- NA
      geneTag <- NA
      
      # Note its type
      type <- cols[2]
      
      # Remove ">" if present in coordinate information
      cols[3] <- gsub(pattern=">", replacement="", cols[3])
      
      # Get its coordinates and note if complement
      if(grepl(cols[3], pattern="complement") == TRUE && 
         grepl(cols[3], pattern="join") == FALSE){
        
        direction <- "COMPLEMENT"
        
        # Remove "complement" and brackets from coordinates
        coordinates <- substr(cols[3], start=12, stop=nchar(cols[3])-1)
        
        # Store the coordinates
        coordinates <- as.numeric(strsplit(coordinates, split="[.]+")[[1]])
      
      # Note that features can be split over fragments
      }else if(grepl(cols[3], pattern="join") == TRUE){
        
        direction <- "COMPLEMENT"
        
        # Feature split over multiple regions
        coordinates <- substr(cols[3], start=17, stop=nchar(cols[3])-2)
        
        # Store the coordinates
        coordinates <- as.numeric(strsplit(coordinates, split="[.]+|,")[[1]])
        
      }else{
        
        # Store the coordinates
        direction <- "FORWARD"
        coordinates <- as.numeric(strsplit(cols[3], split="[.]+")[[1]])
      }
      
    }else if(grepl(cols[2], pattern="/locus_tag") == TRUE){
      
      # Note the locus tag
      locusTag <- gsub(pattern="\"", replacement="",
                       strsplit(cols[2], split="=")[[1]][2])
      locusTag <- toupper(locusTag)

    }else if(grepl(cols[2], pattern="/gene=") == TRUE){
      
      # Note the locus tag
      geneTag <- gsub(pattern="\"", replacement="",
                       strsplit(cols[2], split="=")[[1]][2])
      geneTag <- toupper(geneTag)
      
      # Check whether reached end of feature section  
    }else if(cols[1] == "SQ"){
      foundSequence <- TRUE
      next
    }
    
    ### Sequence ###
    if(foundSequence == TRUE && cols[1] != "//"){
      
      # Get the sequence from current line
      sequenceParts[length(sequenceParts) + 1] <- cols[2]
    }
  }

  # Parse the sequence
  sequence <- paste(sequenceParts, collapse="")
  sequence <- toupper(gsub(pattern=" ", replacement="", sequence))
  
  # Return the feature information and sequence
  output <- list(
    "Features"=features,
    "Sequence"=strsplit(sequence, split="")[[1]]
  )
  
  return(output)
}
