#### Preparation ####

# Load libraries
library(MASS) # glm model with negative binomial distribution
library(basicPlotteR) # alpha colours in plotting

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the working directory
setwd(file.path("~", "Desktop", "HelpingNagwa"))

#### Read in the data ####

# Read in the egg count data
fileName <- "FarmEggCounts_26-11-19.csv"
eggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=TRUE)

# Fix Category names
eggCounts$Category[eggCounts$Category %in% c("Filly", "fily", "Fily")] <- "Filly"
eggCounts$Category[eggCounts$Category %in% c("Gelding", "geloling", "Golding")] <- "Gelding"
eggCounts$Category <- factor(eggCounts$Category)

# Set farm code to factor
eggCounts$Farm.code <- as.factor(eggCounts$Farm.code)

# Note NA values
eggCounts[eggCounts == 999] <- NA
eggCounts$Category <- factor(eggCounts$Category)
eggCounts[eggCounts == ""] <- NA
eggCounts$Gender <- factor(eggCounts$Gender)

#### Plot the data ####

# Plot exploratory plots for the strongyles
generateExploratoryPlots(response="strongyles", 
                         columns=c("Farm.code", "Age", "Category", "Gender", "Species", "parascaris"),
                         data=eggCounts)

# Plot exploratory plots for the strongyles
generateExploratoryPlots(response="parascaris", 
                         columns=c("Farm.code", "Age", "Category", "Gender", "Species", "strongyles"),
                         data=eggCounts)

#### Fit a negative binomial model on the strongyles egg counts ####

nbModel <- glm.nb(strongyles ~ Farm.code + Species + Gender + Age + Category, data=eggCounts)

#### Fit an ANOVA model on the strongyles egg counts ####

anovaModel <- aov(strongyles ~ Farm.code + Gender + Category, data=eggCounts)

#### FUNCTIONS ####

generateExploratoryPlots <- function(response, columns, data){
  
  # Examine each of the columns
  for(column in columns){
    
    # Check if current column is a factor
    if(is.factor(data[, column])){
      
      # Produce a boxplot
      boxplot(eggCounts[, response] ~ eggCounts[, column], las=1, frame=FALSE, outcol=rgb(0,0,0, 0),
              ylab="Faecal egg count (eggs/gram)", main=response, xlab=column)
      spreadPointsMultiple(data=eggCounts, responseColumn=response, categoriesColumn=column, plotOutliers=TRUE,
                           pointCex=0.5, col="red")
      
    }else{
      
      # Plot a scatter plot
      plot(x=eggCounts[, column], y=eggCounts[, response], las=1, pch=19, col=rgb(0,0,0, 0.1),
           ylab="Faecal egg count (eggs/gram)", main=response, xlab=column, bty="n")
    }
  }
}

