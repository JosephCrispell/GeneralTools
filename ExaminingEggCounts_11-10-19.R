#### Preparation ####

# Load libraries
library(basicPlotteR)
library(randomForest)

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the working directory
setwd(file.path("~", "Desktop", "HelpingNagwa"))

#### Read in the data ####

# Read in the corrected egg count data (eggs per gram)
fileName <- "NagwaAnalysis_EPG_11-10-19.csv"
correctedEggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the raw egg count data
fileName <- "NagwaAnalysis_EC_25-10-19.csv"
rawEggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

#### Plot a data summary ####

# Plot the raw data
pdf(paste0("NagwaAnalysis_Figures_", date, ".pdf"))
plotRawEggCounts(correctedEggCounts, yLab="Faecal egg count (eggs/gram)")
plotRawEggCounts(rawEggCounts, yLab="Faecal egg count (raw)")
dev.off()

#### Transform table structure for use in statistical analyses ####

# Transform the table into long form
correctedEggCounts <- transformTableIntoLongform(correctedEggCounts)

#### Fit linear model ####

# Fit a linear model to determine whether the methods are giving significantly results
linearModel <- glm(EggCount ~ Method + Sample + Replicate, family=poisson(), data=correctedEggCounts)

# Summarise the linear model results
summary(linearModel)

#### Fit random forest regression model ###

# Fit a random forest model
rfModel <- randomForest(correctedEggCounts$EggCount ~ ., data=correctedEggCounts[, -1], mtry=3, importance=TRUE, ntree=1000, 
                        keep.forest=TRUE, norm.votes=FALSE, proximity=FALSE, do.trace=FALSE)

# Examine whether number of trees was sufficient (should plateau)
plot(rfModel, las=1)

# Calculate the R squared value
rSq <- rfModel$rsq[length(rfModel$rsq)]

# Look at the variable importance
rfModel$importance

#### Fit one-way ANOVA ####

# Lot's of helpful information here: http://www.sthda.com/english/wiki/one-way-anova-test-in-r

# Fit a one-way anova model
oneWayAnova <- aov(EggCount ~ Method, data=correctedEggCounts)

# Summarise the one-way anova results
summary(oneWayAnova)

#### Fit two-way ANOVA ####

# Lot's of helpful information here: http://www.sthda.com/english/wiki/two-way-anova-test-in-r

# Fit a one-way anova model
twoWayAnova <- aov(EggCount ~ Method + Sample, data=correctedEggCounts)

# Summarise the one-way anova results
summary(twoWayAnova)

#### Fit two-way ANOVA for each sample seperately ####

# Sample 1
anovaSample1 <- aov(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 1, ])
summary(anovaSample1)

# Sample 2
anovaSample2 <- aov(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 2, ])
summary(anovaSample2)

# Sample 3
anovaSample3 <- aov(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 3, ])
summary(anovaSample3)

#### Mann Whitney U test ####

# Fit a pairwise Mann Whitney U test and correct for multiple testing
mannWhitney <- pairwise.wilcox.test(x=correctedEggCounts$EggCount,
                                    g=correctedEggCounts$Method, p.adjust.method="bonferroni")
mannWhitney

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
      output[nrow(output), "Replicate"] <- paste0(sample, ".", replicate)
    }
  }
  
  # Convert the method and sample columns to factors
  output$Method <- as.factor(output$Method)
  output$Sample <- as.factor(output$Sample)
  output$Replicate <- as.factor(output$Replicate)
  
  return(output)
}

plotRawEggCounts <- function(eggCounts, yLab){
  
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
  plot(x=NULL, y=NULL, xlim=c(0.5, 3.5), ylim=range(eggCounts[, -1], na.rm=TRUE), las=1, bty="n", xaxt="n",
       ylab=yLab, xlab="",
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
      points(x=c(xPosition, xPosition), y=range(values, na.rm=TRUE), type="l", col=methodColours[[method]], xpd=TRUE)
      points(x=c(xPosition, xPosition, xPosition), y=values, pch=19, col=setAlpha(methodColours[[method]], 0.5), xpd=TRUE)
    }
  }
  
  # Add a legend
  legend("topright", legend=names(methodColours), text.col=unlist(methodColours), bty="n")
  
  # Reset the plotting margins
  par(mar=currentMar)
}
