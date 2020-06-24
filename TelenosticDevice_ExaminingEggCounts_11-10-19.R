#### Preparation ####

# Load libraries
library(MASS) # glm model with negative binomial distribution
library(basicPlotteR) # alpha colours in plotting
library(broom) # simple table from model outputs

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the working directory
setwd(file.path("~", "Desktop", "HelpingNagwa", "ComparingEggCountsMethods"))

#### Read in the data ####

# Read in the corrected egg count data (eggs per gram)
fileName <- "NagwaAnalysis_EPG_11-10-19.csv"
correctedEggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the raw egg count data
fileName <- "NagwaAnalysis_EC_25-10-19.csv"
rawEggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

#### Plot a data summary ####

# Plot the raw data
jpeg(paste0("NagwaAnalysis_Figures_corrected_", date, "_%03d.jpeg"))
plotRawEggCounts(correctedEggCounts, yLab="Faecal egg count (eggs/gram)", blackAndWhite=TRUE)
plotRawEggCounts(correctedEggCounts, yLab="Faecal egg count (eggs/gram)", blackAndWhite=TRUE)
dev.off()
jpeg(paste0("NagwaAnalysis_Figures_raw_", date, "_%03d.jpeg"))
plotRawEggCounts(rawEggCounts, yLab="Faecal egg count (raw)")
plotRawEggCounts(rawEggCounts, yLab="Faecal egg count (raw)", blackAndWhite=TRUE)
dev.off()

#### Transform table structure for use in statistical analyses ####

# Transform the table into long form
correctedEggCounts <- transformTableIntoLongform(correctedEggCounts)

#### Fit linear model ####

# Fit a linear model to determine whether the methods are giving significantly results
linearModel <- glm.nb(EggCount ~ Method + Sample + Replicate, data=correctedEggCounts)

# Summarise the linear model results
summary(linearModel)

#### Fit the linear model for each sample separately ####

# Sample 1
linearModelSample1 <- glm.nb(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 1, ])
summary(linearModelSample1)
write.table(tidy(linearModelSample1), file=paste("NegativeBinomial_sample1_", date, ".csv"), quote=FALSE, sep=",",
            row.names=FALSE)

# Sample 2
linearModelSample2 <- glm.nb(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 2, ])
summary(linearModelSample2)
write.table(tidy(linearModelSample2), file=paste("NegativeBinomial_sample2_", date, ".csv"), quote=FALSE, sep=",",
            row.names=FALSE)

# Sample 3
linearModelSample3 <- glm.nb(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 3, ])
summary(linearModelSample3)
write.table(tidy(linearModelSample3), file=paste("NegativeBinomial_sample3_", date, ".csv"), quote=FALSE, sep=",",
            row.names=FALSE)

#### Fit two-way ANOVA ####

# Lot's of helpful information here: http://www.sthda.com/english/wiki/two-way-anova-test-in-r
# Assumptions of ANOVA:
# - The observations are obtained independently and randomly from the population defined by the factor levels
# - The data of each factor level are normally distributed.
# - These normal populations have a common variance. (Levene’s test can be used to check this.)

# Fit a two-way anova model
twoWayAnova <- aov(EggCount ~ Method + Sample, data=correctedEggCounts)

# Summarise the one-way anova results
summary(twoWayAnova)

#### Fit two-way ANOVA for each sample seperately ####

# Sample 1
anovaSample1 <- aov(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 1, ])
summary(anovaSample1)
write.table(tidy(anovaSample1), file=paste("ANOVA_sample1_", date, ".csv"), quote=FALSE, sep=",",
            row.names=FALSE)

# Sample 2
anovaSample2 <- aov(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 2, ])
summary(anovaSample2)
write.table(tidy(anovaSample2), file=paste("ANOVA_sample2_", date, ".csv"), quote=FALSE, sep=",",
            row.names=FALSE)

# Sample 3
anovaSample3 <- aov(EggCount ~ Method + Replicate, data=correctedEggCounts[correctedEggCounts$Sample == 3, ])
summary(anovaSample3)
write.table(tidy(anovaSample3), file=paste("ANOVA_sample3_", date, ".csv"), quote=FALSE, sep=",",
            row.names=FALSE)

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

plotRawEggCounts <- function(eggCounts, yLab, blackAndWhite=FALSE){
  
  # Note the padding locations for the analyses
  padding <- seq(from=-0.1, to=0.1, length.out=4)

  # Note the colour for each method
  methodStyle <- list(
    "MC"=list(col="red", pch=19, pointAlpha=0.5, lty=1, 
              name="McMaster", padding=padding[1], pointCex=1),
    "MiniFLOTAC"=list(col="green", pch=19, pointAlpha=0.5, lty=1, 
                      name="MiniFLOTAC", padding=padding[2], pointCex=1), 
    "TelenosticManual"=list(col="blue", pch=19, pointAlpha=0.5, lty=1, 
                            name="TelenosticManual", padding=padding[3], pointCex=1),
    "TelenosticAutomatic"=list(col="black", pch=19, pointAlpha=0.5, lty=1, 
                               name="TelenosticAutomatic", padding=padding[4], pointCex=1)
  )
  
  # Change styles if plot is black and white
  if(blackAndWhite){
    methodStyle <- list(
      "MC"=list(col="black", pch=15, pointAlpha=1, lty=5, 
                name="McMaster", padding=padding[1], pointCex=1),
      "MiniFLOTAC"=list(col="black", pch=17, pointAlpha=1, lty=4, 
                        name="MiniFLOTAC", padding=padding[2], pointCex=1), 
      "TelenosticManual"=list(col="darkgrey", pch=18, pointAlpha=1, lty=3, 
                              name="TelenosticManual", padding=padding[3], pointCex=1.5),
      "TelenosticAutomatic"=list(col="darkgrey", pch=19, pointAlpha=1, lty=1, 
                                 name="TelenosticAutomatic", padding=padding[4], pointCex=1)
    )
  }
 
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
    for(method in names(methodStyle)){
      
      # Note the X location
      xPosition <- sampleXPositions[[eggCounts[row, "Sample"]]] + methodStyle[[method]]$padding
      
      # Get the yValues
      values <- eggCounts[row, grepl(colnames(eggCounts), pattern=method)]
      
      # Plot the data for the current sample
      points(x=c(xPosition, xPosition), y=range(values, na.rm=TRUE), 
             type="l", lty=methodStyle[[method]]$lty,
             col=methodStyle[[method]]$col, xpd=TRUE)
      points(x=c(xPosition, xPosition, xPosition), y=values, 
             pch=methodStyle[[method]]$pch, cex=methodStyle[[method]]$pointCex,
             col=setAlpha(methodStyle[[method]]$col, methodStyle[[method]]$pointAlpha), xpd=TRUE)
    }
  }
  
  # Add a legend
  legend("topright", 
         legend=sapply(names(methodStyle), 
                       FUN=function(method, methodStyle){
                          return(methodStyle[[method]]$name)
                        }, methodStyle), 
         text.col=sapply(names(methodStyle), 
                         FUN=function(method, methodStyle){
                           return(methodStyle[[method]]$col)
                         }, methodStyle),
         bty="n")
  if(blackAndWhite){
    legend("topright", 
           legend=sapply(names(methodStyle), 
                         FUN=function(method, methodStyle){
                           return(methodStyle[[method]]$name)
                         }, methodStyle), 
           text.col=sapply(names(methodStyle), 
                           FUN=function(method, methodStyle){
                             return(methodStyle[[method]]$col)
                           }, methodStyle),
           pch=sapply(names(methodStyle), 
                      FUN=function(method, methodStyle){
                        return(methodStyle[[method]]$pch)
                      }, methodStyle),
           lty=sapply(names(methodStyle), 
                      FUN=function(method, methodStyle){
                        return(methodStyle[[method]]$lty)
                      }, methodStyle),
           col=sapply(names(methodStyle), 
                      FUN=function(method, methodStyle){
                        return(methodStyle[[method]]$col)
                      }, methodStyle),
           bty="n")
  }
  
  # Reset the plotting margins
  par(mar=currentMar)
}
