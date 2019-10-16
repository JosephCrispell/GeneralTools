# Load libraries
library(basicPlotteR)
library(randomForest)

# Read in the data
fileName <- "/home/josephcrispell/Desktop/HelpingNagwa/AnalysisEPGandECs_11-10-19.csv"
eggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Plot the raw data
pdf(paste0(substr(fileName, 1, nchar(fileName)-4), ".pdf"))
plotRawEggCounts(eggCounts)
dev.off()

# Transform the table into long form
eggCounts <- transformTableIntoLongform(eggCounts)

# Fit a linear model to determine whether the methods are giving significantly results
# Sample and Replicate should be a single variable
modelFit <- glm(EggCount ~ Method + Sample + Replicate, family=poisson(), data=eggCounts)
summary(modelFit)

# Fit a random forest model
# Fit the Random Forest model
rfModel <- randomForest(eggCounts$EggCount ~ ., data=eggCounts[, -1], mtry=3, importance=TRUE, ntree=1000, 
                        keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE, do.trace=FALSE)
rSq <- rfModel$rsq[length(rfModel$rsq)]
#plot(rfModel, las=1)
rfModel$importance


#### FUNCTIONS ####

transformTableIntoLongform <- function(eggCounts){
  
  # Initialise a table to store the data in long format
  output <- data.frame("EggCount"=c(), "Method"=c(), "Sample"=c(), "Replicate"=c())
  
  # Examine every row in the table
  for(row in 1:nrow(eggCounts)){
    
    # Get the current sample number
    sample <- ifelse(grepl(eggCounts[row, "Sample"], pattern="Sample 1"), 1, 
                     ifelse(grepl(eggCounts[row, "Sample"], pattern="Sample 2"), 2, 3))
    
    # Get the current replicate number
    replicate <- strsplit(eggCounts[row, "Sample"], split="\\.")[[1]][2]
    
    # Examine every column in the table
    for(col in 2:ncol(eggCounts)){
      
      # Get the current method
      method <- strsplit(colnames(eggCounts)[col], split="_")[[1]][1]
      
      # Store the current value in the table
      output[nrow(output) + 1, "EggCount"] <- eggCounts[row, col]
      output[nrow(output), "Method"] <- method
      output[nrow(output), "Sample"] <- sample
      output[nrow(output), "Replicate"] <- replicate
    }
  }
  
  # Convert the method and sample columns to factors
  output$Method <- as.factor(output$Method)
  output$Sample <- as.factor(output$Sample)
  output$Replicate <- as.factor(output$Replicate)
  
  return(output)
}

plotRawEggCounts <- function(eggCounts){
  
  # Note the padding locations for the analyses
  padding <- seq(from=-0.1, to=0.1, length.out=4)
  methodPadding <- list("MC"=padding[1], "MiniFLOTAC"=padding[2], "TelenosticManual"=padding[3], "TelenosticAutomatic"=padding[4])
  
  # Note the colour for each method
  methodColours <- list("MC"="red", "MiniFLOTAC"="green", "TelenosticManual"="blue", "TelenosticAutomatic"="black")
  
  # Note the sample replicate locations
  sampleXPositions <- list("Sample 1.1"=0.66, "Sample 1.2"=1, "Sample 1.3"=1.33,
                           "Sample 2.1"=1.66, "Sample 2.2"=2, "Sample 2.3"=2.33,
                           "Sample 3.1"=2.66, "Sample 3.2"=3, "Sample 3.3"=3.33)
  
  # Set the plotting margins
  currentMar <- par()$mar
  par(mar=c(9,4.1,4.1,2.1))
  
  # Create an empty plot
  plot(x=NULL, y=NULL, xlim=c(0.5, 3.5), ylim=range(eggCounts[, c(2,3,4)]), las=1, bty="n", xaxt="n",
       ylab="Faecal egg count (eggs/gram)", xlab="",
       main="Comparing methods to estimate faecal egg counts")
  
  # Add an X axis
  axis(side=1, at=c(0.66, 1, 1.33), labels=c("A.1", "A.2", "A.3"))
  axis(side=1, at=c(1.66, 2, 2.33), labels=c("B.1", "B.2", "B.3"))
  axis(side=1, at=c(2.66, 3, 3.33), labels=c("C.1", "C.2", "C.3"))
  mtext("Replicate", side=1, line=2.5)
  axis(side=1, at=c(1,2,3), labels=c("A", "B", "C"), line=5)
  mtext("Sample", side=1, line=7.5)
  
  # Examine every row of the table
  for(row in 1:nrow(eggCounts)){
    
    # Examine each analysis
    for(method in names(methodPadding)){
      
      # Note the X location
      xPosition <- sampleXPositions[[eggCounts[row, "Sample"]]] + methodPadding[[method]]
      
      # Get the yValues
      values <- eggCounts[row, grepl(colnames(eggCounts), pattern=method)]
      
      # Plot the data for the current sample
      points(x=c(xPosition, xPosition), y=range(values), type="l", col=methodColours[[method]])
      points(x=c(xPosition, xPosition, xPosition), y=values, pch=19, col=setAlpha(methodColours[[method]], 0.5))
    }
  }
  
  # Add a legend
  legend("topright", legend=names(methodColours), text.col=unlist(methodColours), bty="n")
  
  # Reset the plotting margins
  par(mar=currentMar)
}
