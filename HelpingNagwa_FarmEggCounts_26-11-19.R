#### Preparation ####

# Load libraries
library(MASS) # glm model with negative binomial distribution
library(basicPlotteR) # alpha colours in plotting
library(broom) # simple table from model outputs

# Get the current date
date <- format(Sys.Date(), "%d-%m-%y")

# Set the working directory
setwd(file.path("~", "Desktop", "HelpingNagwa", "ExamineFarmEggCounts"))

#### Read in the data ####

# Read in the horse egg count data
fileName <- "HorseFarmEggCounts_21-12-19.csv"
horseEggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Read in the donkey egg count data
fileName <- "DonkeyFarmEggCounts_21-12-19.csv"
donkeyEggCounts <- read.table(fileName, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Tidy up both datasets
horseEggCounts <- tidyData(horseEggCounts)
donkeyEggCounts <- tidyData(donkeyEggCounts)

#### Plot the data ####

# Plot the egg count distribution for horses and donkeys
pdf(paste0("EggCountDistributions_", date, ".pdf"), width=14, height=14)
par(mfrow=c(2,2))
hist(horseEggCounts$strongyles, las=1, breaks=20, col="red", main="Strongyles egg counts in horses",
     xlab="Faecal egg count (eggs/gram)", cex.axis=1.5, cex.lab=1.5, cex.main=2)
hist(horseEggCounts$parascaris, las=1, breaks=20, col="red", main="Parascaris egg counts in horses",
     xlab="Faecal egg count (eggs/gram)", cex.axis=1.5, cex.lab=1.5, cex.main=2)
hist(donkeyEggCounts$strongyles, las=1, breaks=20, col="red", main="Strongyles egg counts in donkeys",
     xlab="Faecal egg count (eggs/gram)", cex.axis=1.5, cex.lab=1.5, cex.main=2)
hist(donkeyEggCounts$parascaris, las=1, breaks=20, col="red", main="Strongyles egg counts in donkeys",
     xlab="Faecal egg count (eggs/gram)", cex.axis=1.5, cex.lab=1.5, cex.main=2)
dev.off()

# Plot exploratory plots for the strongyles in horses
pdf(paste0("ExploratoryPlots_HORSES_", date, ".pdf"))
generateExploratoryPlots(response="strongyles", 
                         columns=c("Farm.code", "Category", "TimeSinceTreatment", "Season",
                                   "Type_treatment"),
                         data=horseEggCounts)
dev.off()

# Plot exploratory plots for the strongyles in donkeys
pdf(paste0("ExploratoryPlots_DONKEYS_", date, ".pdf"))
generateExploratoryPlots(response="strongyles", 
                         columns=c("Farm.code", "Category", "TimeSinceTreatment", "Season",
                                   "Type_treatment"),
                         data=donkeyEggCounts)
dev.off()

# Compare the egg counts in donkeys and horses
pdf(paste0("ComparingDonkeysAndHorses_", date, ".pdf"))
generateExploratoryPlots(response="strongyles", columns=c("Species"),
                         data=rbind(horseEggCounts, donkeyEggCounts), changeMargins=FALSE)
generateExploratoryPlots(response="parascaris", columns=c("Species"),
                         data=rbind(horseEggCounts, donkeyEggCounts), changeMargins=FALSE)
dev.off()

#### Fit a negative binomial model on the strongyles egg counts ####

# Fit for the horses data
nbModel <- glm.nb(strongyles ~ Farm.code + Category + TimeSinceTreatment + Season + Type_treatment,
                  data=horseEggCounts, maxit=50)
write.table(tidy(nbModel), file=paste("NegativeBinomial_Strongyles_HORSES_", date, ".csv"),
            quote=FALSE, sep=",", row.names=FALSE)
summary(nbModel)

# Fit for the donkey data
nbModel <- glm.nb(strongyles ~ Farm.code + Category + TimeSinceTreatment + Season,
                  data=horseEggCounts, maxit=50)
write.table(tidy(nbModel), file=paste("NegativeBinomial_Strongyles_DONKEYS_", date, ".csv"),
            quote=FALSE, sep=",", row.names=FALSE)
summary(nbModel)

#### FUNCTIONS ####

tidyData <- function(eggCounts){

  # Edit the drug types
  eggCounts$Type_treatment[eggCounts$Type_treatment == "Pyrantel "] <- "Pyrantel"
  eggCounts$Type_treatment[eggCounts$Type_treatment == "5day Benzimidazolee"] <- "Benzimidazole"
  eggCounts$Type_treatment[eggCounts$Type_treatment == "Moxidectin "] <- "Moxidectin"
  eggCounts$Type_treatment[eggCounts$Type_treatment == "Pancur and Moxidectin "] <- "Pancur and Moxidectin"
    
  # Identify the NA values
  eggCounts[eggCounts == "Missing"] <- NA
  eggCounts[eggCounts == ""] <- NA
  
  # Create the factors
  eggCounts$Farm <- as.factor(eggCounts$Farm)
  eggCounts$Farm.code <- as.factor(eggCounts$Farm.code)
  eggCounts$Age <- as.numeric(eggCounts$Age)
  eggCounts$Category <- as.factor(eggCounts$Category)
  eggCounts$Gender <- as.factor(eggCounts$Gender)
  eggCounts$Species <- as.factor(eggCounts$Species)
  eggCounts$strongyles <- as.numeric(eggCounts$strongyles)
  eggCounts$parascaris <- as.numeric(eggCounts$parascaris)
  eggCounts$Type_treatment <- as.factor(eggCounts$Type_treatment)
  eggCounts$Season <- as.factor(eggCounts$Season)
  
  # Note the dates - 07-Dec-14
  eggCounts$Data_last_treat <- as.Date(eggCounts$Data_last_treat, format="%d-%b-%y")
  eggCounts$Data_test <- as.Date(eggCounts$Data_test, format="%d-%b-%y")
  
  # Calculate the time since last treament
  eggCounts$TimeSinceTreatment <- as.numeric(eggCounts$Data_test - eggCounts$Data_last_treat)
  
  # Check for negative values in time since last treatment
  wrongDateIndices <- which(eggCounts$TimeSinceTreatment < 0)
  if(length(wrongDateIndices)){
    warning("Incorrect dates found in ", length(wrongDateIndices), " rows: \n", 
            paste(wrongDateIndices, collapse=", "), "\n", "replaced with NAs")
    eggCounts$TimeSinceTreatment[wrongDateIndices] <- NA
  }
  
  return(eggCounts)
}

generateExploratoryPlots <- function(response, columns, data, changeMargins=TRUE, ...){

  # Examine each of the columns
  for(column in columns){
    
    # Check if current column is a factor
    if(is.factor(data[, column])){
      
      # Set the plotting margins
      currentMar <- par()$mar
      if(changeMargins){
        par(mar=c(10, 4.1, 4.1, 2.1))
      }
      
      # Produce a boxplot
      boxplot(data[, response] ~ data[, column], las=2, frame=FALSE, outcol=rgb(0,0,0, 0),
              ylab="Faecal egg count (eggs/gram)", main=response, 
              xlab=ifelse(length(unique(data[, column])) == 1, 
                          as.character(unique(data[, column])[1]), ""), ...)
      spreadPointsMultiple(data=data, responseColumn=response, categoriesColumn=column, plotOutliers=TRUE,
                           pointCex=0.5, col="red")
      
      # Reset the plotting margins
      par(mar=currentMar)
    }else{
      
      # Plot a scatter plot
      plot(x=data[, column], y=data[, response], las=1, pch=19, col=rgb(0,0,0, 0.1),
           ylab="Faecal egg count (eggs/gram)", main=response, xlab=column, bty="n", ...)
    }
  }
}

