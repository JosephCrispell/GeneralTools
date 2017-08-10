############
# Set path #
############

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

#####################################
# Read in variant position coverage #
#####################################

# Read in the genome coverage file
coverageFile <- paste(path, "vcfFiles/IsolateVariantPositionCoverage_01-08-2017.txt", sep="")
coverage <- read.table(coverageFile, header=TRUE, stringsAsFactors=FALSE)

# Parse the Isolate column
coverage$Isolate <- parseIsolateColumn(coverage$Isolate)

#############################
# Plot the isolate coverage #
#############################

# Open a pdf
file <- paste(substr(coverageFile, 1, nchar(coverageFile) - 4), ".pdf", sep="")
pdf(file)

plot(coverage$VariantPositionCoverage, pch=20, xaxt="n", xlab="", bty="n", las=1,
     col=ifelse(grepl(x=coverage$Isolate, pattern="WB"), rgb(1,0,0, 0.5),
                rgb(0,0,1, 0.5)),
     ylab="Proportion", main="Variant Position Coverage")

text(y=coverage$VariantPositionCoverage,
     x=1:nrow(coverage),
     labels = coverage$Isolate, cex=0.1,
     col=ifelse(coverage$VariantPositionCoverage < 0.75, rgb(0,0,0, 1), rgb(0,0,0, 0)))

coverage$Species <- "COW"
coverage$Species[grepl(x=coverage$Isolate, pattern="WB") == TRUE] <- "BADGER"
coverage$Species <- as.factor(coverage$Species)

boxplot(coverage$VariantPositionCoverage ~ coverage$Species, 
        ylim=c(0,1), ylab="Proportion", border=c("red", "blue"),
        names=c("Badgers", "Cattle"), outline=FALSE,
        las=1, pch=20, main="Isolate Variant Position Coverage")

stripchart(coverage$VariantPositionCoverage ~ coverage$Species,
           vertical = TRUE, jitter=0.2,
           method = "jitter", add = TRUE, pch = 21,
           col = c(rgb(1,0,0, 0.5), rgb(0,0,1, 0.5)),
           bg=rgb(0.5,0.5,0.5, 0.5))



###################################
# Note the poor coverage isolates #
###################################

# Set coverage threshold
threshold <- 0.75
abline(h=threshold, lty=2, lwd=2)
dev.off()

# Select the isolates with coverage above threshold
selected <- coverage[coverage$VariantPositionCoverage >= threshold, ]
cat(paste(length(which(grepl(x=selected$Isolate, pattern="WB"))), " badgers kept"))
cat(paste(length(which(grepl(x=selected$Isolate, pattern="TB"))), " cattle kept"))

# Select those below
remove <- coverage[coverage$VariantPositionCoverage < threshold, ]
cat(paste(length(which(grepl(x=remove$Isolate, pattern="WB"))), " badgers removed"))
cat(paste(length(which(grepl(x=remove$Isolate, pattern="TB"))), " cattle removed"))

# Note the isolates to remove
isolatesToRemove <- data.frame(Isolate=remove$Isolate, 
                               Reason=rep("PoorVariantPositionCoverage", nrow(remove)),
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

# Examine culture IDs - whether multiple present choose isolate with highest coverage
resequencedIsolatesToRemove <- examineCultureIsolatesOfCattleChooseWhenResequenced(
  cultureIDs, coverage)

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

# Get the date from the variant position coverage file name
parts <- strsplit(strsplit(coverageFile, split=".tx")[[1]][1], split="_")[[1]]
date <- parts[length(parts)]

# Create the output file
file <- paste(path, "vcfFiles/isolatesToRemove_VPCoverage-Resequenced_", date, ".txt",
              sep="")
write.table(isolatesToRemove, file, quote=FALSE, row.names=FALSE, sep="\t")

#############
# FUNCTIONS #
#############

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
                          "VariantPositionCoverage"]
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
