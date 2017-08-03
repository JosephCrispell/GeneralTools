################
# Set the path #
################

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Review/"

############################
# Read in the Scholar Hits #
############################

scholarHitsFile <- paste(path, "Mycobacterium_GenomeSequence_2000.csv", sep="")
scholarHits <- read.table(scholarHitsFile, header=TRUE, stringsAsFactors=FALSE, sep=",")

#################
# Plot the data #
#################

# Define year range of scholar search
yearRange <- c(2000, 2017)

# Get year counts
yearCounts <- countHitsPerYear(scholarHits, yearRange)
plot(y=yearCounts, x=yearRange[1]:yearRange[2], type="o", las=1,
     xlab="Year", ylab= "Number of Publications",
     main="N. WGS articles published per year for MTBC")


#############
# FUNCTIONS #
#############

countHitsPerYear <- function(scholarHits, yearRange){
  
  counts <- rep(0, (yearRange[2] - yearRange[1]) + 1)
  
  for(row in 1:nrow(scholarHits)){
    
    if(scholarHits[row, "Year"] >= yearRange[1]){
      yearIndex <- (scholarHits[row, "Year"] - yearRange[1]) + 1
      counts[yearIndex] <- counts[yearIndex] + 1
    }
  }
  
  return(counts)
}