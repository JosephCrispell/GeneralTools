##################
# Load Libraries #
##################

library(plotrix) # Draw circle

########
# Path #
########

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

#####################################
# Get a list of isolates grown well #
#####################################

# Read in the table
file <- paste(path, "NewCattle/",
              "Rowlands second list April 2017v 16.11.17.csv", sep="")
isolatesGrown <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

##################################################
# Get spoligotype information for these isolates #
##################################################

# Read in the spoligotype table
file <- paste(path, "Spoligo_APHA_2013-09-09.csv", sep="")
spoligotypeInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

###################
# Make some plots #
###################

# Open a PDF
file <- paste(path, "NewCattle/Sampling_17-11-17.pdf", sep="")
pdf(file)

# Set the margins
par(mar=c(5.1, 4.1, 4.1, 2.1))

subset <- "ALL"
spoligotypeInfoSubset <- selectSubset(isolatesGrown, spoligotypeInfo, subset)
title <- paste(subset, " (n = ", nrow(spoligotypeInfoSubset), ")", sep="")
barplot(table(spoligotypeInfoSubset$Year), las=2, ylab="N. samples", 
        main=title)
plotSpatialRange(title, spoligotypeInfoSubset)

subset <- "GROWN"
spoligotypeInfoSubset <- selectSubset(isolatesGrown, spoligotypeInfo, subset)
title <- paste(subset, " (n = ", nrow(spoligotypeInfoSubset), ")", sep="")
barplot(table(spoligotypeInfoSubset$Year), las=2, ylab="N. samples", 
        main=title)
plotSpatialRange(title, spoligotypeInfoSubset)

subset <- "NOT GROWN"
spoligotypeInfoSubset <- selectSubset(isolatesGrown, spoligotypeInfo, subset)
title <- paste(subset, " (n = ", nrow(spoligotypeInfoSubset), ")", sep="")
barplot(table(spoligotypeInfoSubset$Year), las=2, ylab="N. samples", 
        main=title)
plotSpatialRange(title, spoligotypeInfoSubset)

dev.off()

#############
# FUNCTIONS #
#############

selectSubset <- function(isolatesGrown, spoligotypeInfo, subset){

  # Select subset of isolates
  isolatesGrownSubset <- isolatesGrown
  if(subset == "GROWN"){
    isolatesGrownSubset <- isolatesGrown[
      isolatesGrown$Isolate.grown.16.11.17 == "yes", ]
  }else if(subset == "NOT GROWN"){
    isolatesGrownSubset <- isolatesGrown[
      isolatesGrown$Isolate.grown.16.11.17 != "yes", ]
  }else if(subset != "ALL"){
    cat(paste("Subset selection not recognised:", subset, "\n"))
  }
  
  # Get the spoligotype info for the grown isolate
  spoligotypeInfoSubset <- spoligotypeInfo[
    spoligotypeInfo$Sample.Ref %in% isolatesGrownSubset$Sample.Ref, ]
  
  return(spoligotypeInfoSubset)
}

plotSpatialRange <- function(title, spoligotypeInfoSubset){
  
  rangeX <- range(spoligotypeInfoSubset$Mapx)
  rangeY <- range(spoligotypeInfoSubset$Mapy)
  
  plot(spoligotypeInfoSubset$Mapx, spoligotypeInfoSubset$Mapy, pch=20, col=rgb(0,0,0, 0.5),
       xlim=rangeX, ylim=rangeY,
       xlab=paste(round((rangeX[2] - rangeX[1]) / 1000, digits=2), "km", sep=""),
       ylab=paste(round((rangeY[2] - rangeY[1]) / 1000, digits=2), "km", sep=""),
       xaxt="n", yaxt="n", bty="n", asp=1,
       main=title)
  
  # Add badger centre
  badgerCentre <- c(381761.7, 200964.3)
  points(x=badgerCentre[1], y=badgerCentre[2], col="red", pch=17)
  
  # Add circle around this
  radius <- 10000
  draw.circle(x=badgerCentre[1], y=badgerCentre[2], 
              radius=radius,
              border="black", lty=2)
  text(x=badgerCentre[1], y=badgerCentre[2]+radius,
       labels=paste(round(radius / 1000, digits=2), "km", sep=""),
       pos=3)
}