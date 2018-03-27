############
# Set path #
############

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/"

###############################################################
# Read in variant position coverage - before and after rescue #
###############################################################

# Read in the genome coverage file
file <- paste(path, "vcfFiles/IsolateVariantPositionCoverage_FILTERED_22-03-2018.txt", sep="")
before <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

file <- paste(path, "vcfFiles/IsolateVariantPositionCoverage_RESCUED_22-03-2018.txt", sep="")
after <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Parse the Isolate columns
before$Isolate <- parseIsolateColumn(before$Isolate)
after$Isolate <- parseIsolateColumn(after$Isolate)

#############################
# Plot the isolate coverage #
#############################

# Add species column
after$Species <- "COW"
after$Species[grepl(x=after$Isolate, pattern="WB") == TRUE] <- "BADGER"

# Open a pdf
file <- paste(substr(file, 1, nchar(file) - 4), ".pdf", sep="")
pdf(file)

produceSummaryPlots(before, after)

dev.off()

###################################
# Note the poor coverage isolates #
###################################

# Set the coverage threshold
threshold <- 0.95

# Select the isolates with coverage above threshold
selected <- after[after$Coverage >= threshold, ]
cat(paste(length(which(selected$Species == "BADGER")), " badgers kept"))
cat(paste(length(which(selected$Species == "COW")), " cattle kept"))

# Select those below
remove <- after[after$Coverage < threshold, ]
cat(paste(length(which(remove$Species == "BADGER")), " badgers removed"))
cat(paste(length(which(remove$Species == "COW")), " cattle removed"))

# Note the isolates to remove
isolatesToRemove <- data.frame(Isolate=remove$Isolate, 
                               Reason=rep("Poor variant position coverage following rescue", nrow(remove)),
                               stringsAsFactors=FALSE)

#################################
# Note the resequenced isolates #
#################################

#########
# Badgers

# Note whether WB_id are duplicated
resequencedIsolatesToRemove <- lookForDuplicatedWBsToIdentifyResequencedBadgers(selected)

if(length(resequencedIsolatesToRemove) > 0){
  isolatesToRemove <- addIsolatesToRemove(resequencedIsolatesToRemove, "Resequenced",
                                          isolatesToRemove)
  cat(paste(length(resequencedIsolatesToRemove), " resequenced badgers to remove\n",
            sep=""))
}else{
  cat("No resequenced badger isolates to remove!\n")
}

########
# Cattle

# Read in link table - links sequence ID to original culture
file <- paste(path, "IsolateData/", 
              "LinkTable_StrainID-SequenceNo_20-06-16.csv", sep="")
linkTable <- read.table(file, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Remove sequence IDs that aren't present
linkTable <- linkTable[linkTable$Seq.number %in% selected$Isolate, ]

# Create hashtable for culture IDs
cultureIDs <- createHashtableOfCultureIDs(linkTable)

# Examine culture IDs - whether multiple present choose isolate with highest after
resequencedIsolatesToRemove <- examineCultureIsolatesOfCattleChooseWhenResequenced(
  cultureIDs, after)

if(length(resequencedIsolatesToRemove) > 0){
  isolatesToRemove <- addIsolatesToRemove(resequencedIsolatesToRemove, "Resequenced",
                                          isolatesToRemove)
  cat(paste(length(resequencedIsolatesToRemove), " resequenced cattle to remove\n",
            sep=""))
}else{
  cat("No resequenced cattle isolates to remove!\n")
}

################################################################################
# Create output file noting isolates to ignore - poor coverage and resequenced #
################################################################################

# Create the output file
file <- paste(path, "vcfFiles/isolatesToRemove_22-03-18.txt", sep="")
write.table(isolatesToRemove, file, quote=FALSE, row.names=FALSE, sep="\t")

#############
# FUNCTIONS #
#############

produceSummaryPlots <- function(before, after){
  
  table <- data.frame(Isolate=before$Isolate, stringsAsFactors=FALSE)
  table$Before <- before$Coverage
  table$After <- after$Coverage
  
  values <- c(table$Before, table$After)
  
  plot(table$Before, table$After, las=1, xlab="Before", ylab="After",
       main="Difference between before and after rescue",
       ylim=range(values), xlim=range(values),
       pch=ifelse(grepl(pattern="WB", x=table$Isolate), 19, 17),
       col=ifelse(grepl(pattern="WB", x=table$Isolate), 
                  rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)))
  
  lines(x=range(values), y=range(values), lty=2, col="black")
  for(row in 1:nrow(table)){
    
    lines(x=c(table[row, "Before"], table[row, "Before"]),
          y=c(table[row, "Before"], table[row, "After"]),
          lty=3, lwd=0.5,
          col=ifelse(grepl(pattern="WB", x=table[row, "Isolate"]), 
                     rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)))
  }
  
  legend("bottomright", legend=c("Cow", "Badger"),
         pch=c(17, 19), cex=1, col=c("blue", "red"), 
         text.col=c("blue", "red"), bty='n')
  
  after$Species <- "COW"
  after$Species[grepl(x=after$Isolate, pattern="WB") == TRUE] <- "BADGER"
  after$Species <- as.factor(after$Species)
  
  boxplot(after$Coverage ~ after$Species, 
          ylim=c(0,1), ylab="Proportion", border=c("red", "blue"),
          names=c("Badgers", "Cattle"), outline=FALSE,
          las=1, pch=20, main="Isolate Variant Position Coverage")
  
  stripchart(after$Coverage ~ after$Species,
             vertical = TRUE, jitter=0.2,
             method = "jitter", add = TRUE, pch = 21,
             col = c(rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)),
             bg=rgb(0.5,0.5,0.5, 0.5))
  
#  plot(table$Before, table$After, las=1, xlab="Before", ylab="After",
#       main="Difference between before and after rescue",
#       ylim=range(values), xlim=range(values),
#       pch=ifelse(grepl(pattern="WB", x=table$Isolate), 19, 17),
#       col=ifelse(grepl(pattern="WB", x=table$Isolate), 
#                  rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)))
#  lines(x=range(values), y=range(values), lty=2, col="black")
#  text(table$Before, table$After,
#       labels=table$Isolate, cex=0.5)
}

lookForDuplicatedWBsToIdentifyResequencedBadgers <- function(coverage){
  
  wbIDs <- list()
  
  # Note the sequence IDs for each WB ID
  for(row in 1:nrow(coverage)){
    
    if(grepl(x=coverage[row, "Isolate"], pattern="WB") == TRUE){
      
      id <- strsplit(coverage[row, "Isolate"], split="_")[[1]][1]
      
      if(is.null(wbIDs[[id]]) == TRUE){
        wbIDs[[id]] <- c(coverage[row, "Isolate"])
      }else{
        wbIDs[[id]] <- c(wbIDs[[id]], coverage[row, "Isolate"])
      }
    }
  }
  
  # Examine resequenced WB IDs
  resequencedIsolatesToRemove <- c()
  
  for(wbID in names(wbIDs)){
    
    if(length(wbIDs[[wbID]]) > 1){
      
      # Get the coverage of each isolate
      coverageValues <- getCoverage(wbIDs[[wbID]], coverage)
      
      # Find the index of the max
      maxIndex <- which.max(coverageValues)
      
      # Add the other(s) into array to remove
      for(i in 1:length(wbIDs[[wbID]])){
        if(i != maxIndex){
          resequencedIsolatesToRemove[length(resequencedIsolatesToRemove) + 1] <-
            wbIDs[[wbID]][i]
        }
      }
    }
  }
  
  return(resequencedIsolatesToRemove)
}

examineCultureIsolatesOfCattleChooseWhenResequenced <- function(cultureIDs, coverage){
  
  resequencedIsolatesToRemove <- c()
  
  for(cultureID in names(cultureIDs)){
    
    if(length(cultureIDs[[cultureID]]) > 1){
      
      # Get the coverage of each isolate
      coverageValues <- getCoverage(cultureIDs[[cultureID]], coverage)
      
      # Find the index of the max
      maxIndex <- which.max(coverageValues)
      
      # Add the other(s) into array to remove
      for(i in 1:length(cultureIDs[[cultureID]])){
        if(i != maxIndex){
          resequencedIsolatesToRemove[length(resequencedIsolatesToRemove) + 1] <-
            cultureIDs[[cultureID]][i]
        }
      }
    }
  }
  
  return(resequencedIsolatesToRemove)
  
}

getCoverage <- function(isolates, coverage){
  
  values <- c()
  for(i in 1:length(isolates)){
    
    values[i] <- coverage[which(coverage$Isolate == isolates[i]), 
                          "Coverage"]
  }
  
  return(values)
}

createHashtableOfCultureIDs <- function(linkTable){
  
  cultureIDs <- list()
  for(row in 1:nrow(linkTable)){
    
    if(is.null(cultureIDs[[linkTable[row, "Strain.ID"]]]) == TRUE){
      
      cultureIDs[[linkTable[row, "Strain.ID"]]] <- c(linkTable[row, "Seq.number"])
    }else{
      cultureIDs[[linkTable[row, "Strain.ID"]]] <- c(cultureIDs[[linkTable[row, "Strain.ID"]]],
                                                     linkTable[row, "Seq.number"])
    }
  }
  
  return(cultureIDs)
}

addIsolatesToRemove <- function(ids, reason, table){
  
  for(i in 1:length(ids)){
    
    table[1 + nrow(table), ] <- c(ids[i], reason)
  }
  return(table)
}

parseIsolateColumn <- function(column){
  
  ids <- c()
  for(i in 1:length(column)){
    parts <- strsplit(column[i], split="_")[[1]]
    
    if(grepl(x=column[i], pattern="WB") == TRUE){
      ids[i] <- paste(parts[1], "_", parts[2], sep="")
    }else{
      ids[i] <- parts[1]
    }
    
  }
  
  return(ids)
}
