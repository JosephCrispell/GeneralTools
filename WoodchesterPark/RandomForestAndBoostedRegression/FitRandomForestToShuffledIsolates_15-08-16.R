suppressMessages(library(randomForest))

######################
# Find the Directory #
######################

# Get the path to the necessary files
path <- "/home/josephcrispell/Desktop/Research/Woodchester_CattleAndBadgers/NewAnalyses_22-03-18/Mislabelling/Badger-RF-BR/ShufflingProportion/"

# Set up the file name
prefix <- paste(path, "geneticVsEpiTable_SHUFFLED_", sep="")
date <- "26-03-2018"

# Initialise an array of shuffling proportions
shufflingProps <- c("0.0", "0.05", "0.1", "0.15", "0.2", "0.25", "0.3", "0.35", "0.4",
                    "0.45", "0.5", "0.55", "0.6", "0.65", "0.7", "0.75", "0.8", "0.85",
                    "0.9", "0.95", "1.0")

# How many replicates?
repeats <- seq(from=0, to=9, by=1)

#####################################
# Run a Random Forest On Each Table #
#####################################

# Initialise an output table
nRow <- length(shufflingProps) * length(repeats)
shufflingResults <- data.frame(ProportionShuffled=rep(NA, nRow),
                               Repeat=rep(NA, nRow),
                               PseudoRSquared=rep(NA, nRow),
                               MSE=rep(NA, nRow),
                               Correlation=rep(NA, nRow))
row <- 0
count <- 0
for(prop in shufflingProps){
  count <- count + 1
  
  cat(paste("Examining Sampling Proportion: ", prop, 
            " (", count, "/", length(shufflingProps), ")\t", sep=""))
  
  for(rep in repeats){
    
    # Increment the row
    row <- row + 1
    
    # Build the file name
    file <- paste(prefix, "Prop-", prop, "_", rep, "_", date, ".txt", sep="")
    
    # Read in the table
    table <- read.table(file, header=TRUE)

    # Select only small distances
    table <- table[table$GeneticDistance < 15, ]
    
    # Note the columns to select
    cols <- 2:26
    
    # Fit a Random Forest Model
    nTrees <- 500
    mTry <- 19
    infoRF <- randomForest(table$GeneticDistance~., data=table[, cols],
                           proximity=FALSE, mtry=mTry, importance=FALSE,
                           ntree=nTrees, do.trace=FALSE, keep.forest=FALSE,
                           norm.votes=FALSE)
    
    # Summarise the model fitting
    explained <- infoRF$rsq[length(infoRF$rsq)]
    MSE <- infoRF$mse[length(infoRF$mse)]
    corr <- round(cor(table$GeneticDistance, infoRF$predicted), digits=2)
    
    # Store the results
    shufflingResults[row, 1] <- prop
    shufflingResults[row, 2] <- rep
    shufflingResults[row, 3] <- explained
    shufflingResults[row, 4] <- MSE
    shufflingResults[row, 5] <- corr
    
    cat(".")
  }
  
  cat("\n")
}

shufflingResults$ProportionShuffled <- as.numeric(shufflingResults$ProportionShuffled)
shufflingResults$PseudoRSquared[shufflingResults$PseudoRSquared < 0] <- 0

###############################################
# Read in the data from file - IF ALREADY RUN #
###############################################

file <- paste(path, "ShufflingProportion_RFModelResults_26-03-18.csv", sep="")
shufflingResults <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

####################
# Plot the Results #
####################

# Summarise the results
propValues <- list()
for(row in 1:nrow(shufflingResults)){
  
  if(is.null(propValues[[as.character(shufflingResults[row, 1])]])){
    propValues[[as.character(shufflingResults[row, 1])]] <- rep(NA, 10)
  }
  
  propValues[[as.character(shufflingResults[row, 1])]][shufflingResults[row, 2] + 1] <- shufflingResults[row, 3]
  
}

summary <- data.frame(ProportionShuffled=as.numeric(shufflingProps), 
                      Mean=rep(0, length(shufflingProps)),
                      Median=rep(0, length(shufflingProps)),
                      Min=rep(0, length(shufflingProps)),
                      Max=rep(0, length(shufflingProps)))

for(i in 1:length(shufflingProps)){
  
  values <- propValues[[as.character(as.numeric(shufflingProps[i]))]]

  summary[i, 2] <- mean(values)
  summary[i, 3] <- median(values)
  summary[i, 4] <- min(values)
  summary[i, 5] <- max(values)
}


# Open an output PDF file
file <- paste(path, "ShufflingProportion_26-03-18.pdf", sep="")
pdf(file)

plot(x=NULL, y=NULL,
     xlab="Proportion of Isolates Shuffled",
     ylab="Variation Explained (Psuedo R Squared)",
     main="Effect of Shuffling on a Fitted RF Model",
     las=1, ylim=c(0,1), xlim=c(0,1), bty="n")

for(row in 1:nrow(summary)){
  points(x=c(summary[row, 1], summary[row, 1]),
         y=c(summary[row, 4], summary[row, 5]), 
         type="l", col=rgb(0,0,0, 0.6), lwd=2)
}

points(summary$ProportionShuffled, summary$Mean, type="o", pch=20, cex=1.5)

dev.off()

#############################################
# Print the Fitted RF Model Results to File #
#############################################

file <- paste(path, "ShufflingProportion_RFModelResults_26-03-18.csv", sep="")
write.table(shufflingResults, file, quote=FALSE, row.names=FALSE, sep=",")



