###########################
# Define general settings #
###########################

direction <- "FORWARDS" # FORWARDS or BACKWARDS
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/BASTA/"
factor <- 0.5
cex <- 2

####################################
# Note the dates in the file names #
####################################

deme3SingleOuter_date <- "26-05-17"
deme3_date <- "07-02-17"
deme4_date <- "07-02-17"
deme6_date <- "12-06-17"
deme8_date <- "03-03-17"

#################
# Create figure #
#################

# Open a PDF file
file <- paste(path, "MigrationRateEstimates_21-07-17.pdf", sep="")
pdf(file, height=35, width=14)

# Define plotting grid
par(mfrow=c(5,2))

# Define the line style
par(lend=1) # 0 = round, 1 = butt, 2 = square

# Remove margins
par(mar=c(0,0,0,0))

# 3 Deme
plot3Deme(path=path, popType="Equal", cex=cex,
          date=deme3_date, factor=factor, title="Equal Population Sizes")
plot3Deme(path=path, popType="Varying",cex=cex,
          date=deme3_date, factor=factor, title="Varying Population Sizes")

# 3 Deme - Single Outer
plot3DemeSingleOuter(path=path, popType="Equal", cex=cex,
          date=deme3SingleOuter_date, factor=factor)
plot3DemeSingleOuter(path=path, popType="Varying",cex=cex,
          date=deme3SingleOuter_date, factor=factor)

# 4 Deme
plot4Deme(path=path, popType="Equal", cex=cex,
          date=deme4_date, factor=factor)
plot4Deme(path=path, popType="Varying", cex=cex,
          date=deme4_date, factor=factor)

# 6 Deme
plot6Deme(path=path, popType="Equal", cex=cex,
          date=deme6_date, factor=factor)
plot6Deme(path=path, popType="Varying", cex=cex,
          date=deme6_date, factor=factor)

# 8 Deme
plot8Deme(path=path, popType="Equal", cex=cex,
          date=deme8_date, factor=factor)
plot8Deme(path=path, popType="Varying", cex=cex,
          date=deme8_date, factor=factor)

dev.off()

#############
# FUNCTIONS #
#############

plot3DemeSingleOuter <- function(path, popType, date, factor,
                      cex){
  
  folderName <- paste("3Demes_NoOuterBadgers_", popType, "PopSizes_", date, "/", sep="")
  
  file <- paste(path, folderName,
                "3Demes_NoOuterBadgers_", popType, "PopSizes_", date, ".log", sep="")
  
  logTable <- read.table(file, header=TRUE, sep="\t")
  
  # Remove the burn-in
  burnIn <- round(0.1 * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]
  
  
  #############################################
  # Calculate migration rates forward in time #
  #############################################
  
  logTable <- calculateForwardMigrationRates(logTable)
  
  ###################
  # Migration Rates #
  ###################
  
  # Define line weights
  lineWeights_Rates <- c(
    mean(logTable[, "migModel.forwardRateMatrix_0_1"]),    # 0 -> 1    innerBadger -> innerCow
    mean(logTable[, "migModel.forwardRateMatrix_0_2"]),    # 0 -> 2    innerBadger -> outer
    mean(logTable[, "migModel.forwardRateMatrix_1_0"]),    # 1 -> 0    innerCow -> innerBadger
    mean(logTable[, "migModel.forwardRateMatrix_1_2"]),    # 1 -> 2    innerCow -> outer
    mean(logTable[, "migModel.forwardRateMatrix_2_0"]),    # 2 -> 0    outer -> innerBadger
    mean(logTable[, "migModel.forwardRateMatrix_2_1"])     # 2 -> 1    outer -> innerCow
  )
  
  plot(x=NULL, y=NULL, xlim=c(1, 10), ylim=c(4, 7),
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  
  arrows(x0=7.5, y0=5.25, x1=3.5, y1=5.25, lwd=lineWeights_Rates[3] * factor)    # 1 -> 0    innerCow -> innerBadger
  arrows(x0=8.5, y0=5.5, x1=6.5, y1=6.5, lwd=lineWeights_Rates[4] * factor)        # 1 -> 2    innerCow -> outer
  
  arrows(x0=3.5, y0=4.75, x1=7.5, y1=4.75, lwd=lineWeights_Rates[1] * factor)    # 0 -> 1    innerBadger -> innerCow
  arrows(x0=2, y0=5.5, x1=4, y1=6.5, lwd=lineWeights_Rates[2] * factor)            # 0 -> 2    innerBadger -> outer
  
  arrows(x0=5, y0=6.5, x1=3, y1=5.5, lwd=lineWeights_Rates[5] * factor)            # 2 -> 0    outer -> innerBadger
  arrows(x0=5.5, y0=6.5, x1=7.5, y1=5.5, lwd=lineWeights_Rates[6] * factor)        # 2 -> 1    outer -> innerCow
  
  text(x=c(1.5, 4.75, 7.5), y=c(5, 6.75, 5), 
       labels=c("Inner Badgers", "Outer", "Inner Cattle"), pos=4,
       col=c("red", "black", "blue"), cex=cex)
  
}

plot4Deme <- function(path, popType, date, factor, cex){
  
  # Read in the log file
  folderName <- paste("4Demes_", popType, "PopSizes_", date, "/",
                      sep="")
  
  file <- paste(path, folderName,"4Demes_", popType, "PopSizes_",
                date, ".log", sep="")
  
  logTable <- read.table(file, header=TRUE, sep="\t")
  
  # Remove the burn-in
  burnIn <- round(0.1 * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]
  
  #############################################
  # Calculate migration rates forward in time #
  #############################################
  
  logTable <- calculateForwardMigrationRates(logTable)
  
  ###################
  # Migration Rates #
  ###################
  
  # Define line weights
  lineWeights_Rates <- c(
    mean(logTable[, "migModel.forwardRateMatrix_0_1"]),    # 0 -> 1    innerBadger -> innerCow
    mean(logTable[, "migModel.forwardRateMatrix_0_3"]),    # 0 -> 3    innerBadger -> outerBadger
    mean(logTable[, "migModel.forwardRateMatrix_1_0"]),    # 1 -> 0    innerCow -> innerBadger
    mean(logTable[, "migModel.forwardRateMatrix_1_2"]),    # 1 -> 2    innerCow -> outerCow
    mean(logTable[, "migModel.forwardRateMatrix_2_1"]),    # 2 -> 1    outerCow -> innerCow
    mean(logTable[, "migModel.forwardRateMatrix_2_3"]),    # 2 -> 3    outerCow -> outerBadger
    mean(logTable[, "migModel.forwardRateMatrix_3_0"]),    # 3 -> 0    outerBadger -> innerBadger
    mean(logTable[, "migModel.forwardRateMatrix_3_2"])     # 3 -> 2    outerBadger -> outerCow
  )
  
  plot(x=NULL, y=NULL, xlim=c(1, 10), ylim=c(2,8.5),
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  
  arrows(x0=2, y0=4, x1=2, y1=7, lwd=lineWeights_Rates[1] * factor) # innerBadger -> innerCow
  arrows(x0=3, y0=7, x1=3, y1=4, lwd=lineWeights_Rates[3] * factor) # innerCow -> innerBadger
  
  arrows(x0=8, y0=4, x1=8, y1=7, lwd=lineWeights_Rates[8] * factor) # outerBadger -> outerCow
  arrows(x0=9, y0=7, x1=9, y1=4, lwd=lineWeights_Rates[6] * factor) # outerCow -> outerBadger
  
  arrows(x0=4, y0=3, x1=7.5, y1=3, lwd=lineWeights_Rates[2] * factor) # innerBadger -> outerBadger
  arrows(x0=7.5, y0=4, x1=4, y1=4, lwd=lineWeights_Rates[7] * factor) # outerBadger -> innerBadger
  
  arrows(x0=4, y0=7, x1=7.5, y1=7, lwd=lineWeights_Rates[4] * factor) # innerCow -> outerCow
  arrows(x0=7.5, y0=8, x1=4, y1=8, lwd=lineWeights_Rates[5] * factor) # outerCow -> innerCow
  
  text(x=c(1.5, 1.5, 7.5, 7.5), y=c(3.5, 7.5, 3.5, 7.5), 
       labels=c("Inner Badgers", "Inner Cattle", "Outer Badgers",
                "Outer Cattle"), pos=4, 
       col=c("red", "blue", rgb(1,0,0, 0.5), "blue"), cex=cex)
}

plot3Deme <- function(path, popType, date, factor, title,
                      cex){
  
  folderName <- paste("3Demes_", popType, "PopSizes_", date, "/", sep="")
  
  file <- paste(path, folderName,
                "3Demes_", popType, "PopSizes_", date, ".log", sep="")
  
  logTable <- read.table(file, header=TRUE, sep="\t")
  
  # Remove the burn-in
  burnIn <- round(0.1 * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]

  #############################################
  # Calculate migration rates forward in time #
  #############################################
  
  logTable <- calculateForwardMigrationRates(logTable)
  
  ###################
  # Migration Rates #
  ###################
  
  # Define line weights
  lineWeights_Rates <- c(
    mean(logTable[, "migModel.forwardRateMatrix_0_1"]),    # 0 -> 1    cow -> innerBadger
    mean(logTable[, "migModel.forwardRateMatrix_0_2"]),    # 0 -> 2    cow -> outerBadger
    mean(logTable[, "migModel.forwardRateMatrix_1_0"]),    # 1 -> 0    innerBadger -> cow
    mean(logTable[, "migModel.forwardRateMatrix_1_2"]),    # 1 -> 2    innerBadger -> outerBadger
    mean(logTable[, "migModel.forwardRateMatrix_2_0"]),    # 2 -> 0    outerBadger -> cow
    mean(logTable[, "migModel.forwardRateMatrix_2_1"])     # 2 -> 1    outerBadger -> innerBadger
  )
  
  plot(x=NULL, y=NULL, xlim=c(1, 10), ylim=c(1.5, 10),
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  
  arrows(x0=2, y0=4, x1=4, y1=7, lwd=lineWeights_Rates[3] * factor)         # innerBadger -> cow
  arrows(x0=3.5, y0=2.5, x1=7.5, y1=2.5, lwd=lineWeights_Rates[4] * factor) # innerBadger -> outerBadger
  
  arrows(x0=5, y0=7, x1=3, y1=4, lwd=lineWeights_Rates[1] * factor)         # cow -> innerBadger
  arrows(x0=5.5, y0=7, x1=7.5, y1=4, lwd=lineWeights_Rates[2] * factor)     # cow -> outerBadger
  
  arrows(x0=7.5, y0=3.5, x1=3.5, y1=3.5, lwd=lineWeights_Rates[6] * factor) # outerBadger -> innerBadger
  arrows(x0=8.5, y0=4, x1=6.5, y1=7, lwd=lineWeights_Rates[5] * factor)     # outerBadger -> cow
  
  text(x=c(1.5, 4.75, 7.5), y=c(3, 7.5, 3), 
       labels=c("Inner Badgers", "Cattle", "Outer Badgers"),
       pos=4, col=c("red", "blue", rgb(1,0,0, 0.5)), cex=cex)
  
  text(x=5, y=9, labels=title, cex=cex * 2)
}

plot6Deme <- function(path, popType, date, factor, cex){
  
  folderName <- paste("6Demes_", popType, "PopSizes_", date, "/", sep="")
  
  file <- paste(path, folderName,
                "6Demes_", popType, "PopSizes_", date, ".log", sep="")
  
  logTable <- read.table(file, header=TRUE, sep="\t")
  
  # Remove the burn-in
  burnIn <- round(0.1 * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]

  #############################################
  # Calculate migration rates forward in time #
  #############################################
  
  logTable <- calculateForwardMigrationRates(logTable)
  
  ###################
  # Migration Rates #
  ###################
  
  # BadgerInner       0     SAMPLED
  # CowInner          1     SAMPLED
  # CowOuterEast      2     SAMPLED
  # CowOuterWest      3     SAMPLED
  # BadgerOuterEast   4     UNSAMPLED
  # BadgerOuterWest   5     UNSAMPLED

    # Define line weights
  lineWeights_Rates <- c(
    mean(logTable[, "migModel.forwardRateMatrix_0_1"]), # 1  BadgerInner      -> CowInner
    mean(logTable[, "migModel.forwardRateMatrix_0_4"]), # 2  BadgerInner      -> BadgerOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_0_5"]), # 3  BadgerInner      -> BadgerOuterWest
    
    mean(logTable[, "migModel.forwardRateMatrix_1_0"]), # 4  CowInner         -> BadgerInner
    mean(logTable[, "migModel.forwardRateMatrix_1_2"]), # 5  CowInner         -> CowOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_1_3"]), # 6  CowInner         -> CowOuterWest
    
    mean(logTable[, "migModel.forwardRateMatrix_2_1"]), # 7  CowOuterEast     -> CowInner
    mean(logTable[, "migModel.forwardRateMatrix_2_3"]), # 8  CowOuterEast     -> CowOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_2_4"]), # 9  CowOuterEast     -> BadgerOuterEast
    
    mean(logTable[, "migModel.forwardRateMatrix_3_1"]), # 10 CowOuterWest     -> CowInner
    mean(logTable[, "migModel.forwardRateMatrix_3_2"]), # 11 CowOuterWest     -> CowOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_3_5"]), # 12 CowOuterWest     -> BadgerOuterWest
    
    mean(logTable[, "migModel.forwardRateMatrix_4_0"]), # 13 BadgerOuterEast  -> BadgerInner
    mean(logTable[, "migModel.forwardRateMatrix_4_2"]), # 14 BadgerOuterEast  -> CowOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_4_5"]), # 15 BadgerOuterEast  -> BadgerOuterWest
    
    mean(logTable[, "migModel.forwardRateMatrix_5_0"]), # 16 BadgerOuterWest  -> BadgerInner
    mean(logTable[, "migModel.forwardRateMatrix_5_3"]), # 17 BadgerOuterWest  -> CowOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_5_4"])  # 18 BadgerOuterWest  -> BadgerOuterEast
  )
  
  ### Migration Rates
  plot(x=NULL, y=NULL, xlim=c(0,10), ylim=c(0,10),
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  
  arrows(x0=4.5, y0=4.5, x1=4.5, y1=5.5,
         lwd=lineWeights_Rates[1] * factor)    # 1  BadgerInner      -> CowInner
  arrows(x0=6.5, y0=3.5, x1=8, y1=2.5,
         lwd=lineWeights_Rates[3] * factor)    # 2  BadgerInner      -> BadgerOuterEast
  arrows(x0=4, y0=3.5, x1=2.5, y1=2.5,
         lwd=lineWeights_Rates[2] * factor)    # 3  BadgerInner      -> BadgerOuterWest
  
  arrows(x0=5.5, y0=5.5, x1=5.5, y1=4.5,
         lwd=lineWeights_Rates[4] * factor)    # 4  CowInner         -> BadgerInner
  arrows(x0=6.5, y0=6.5, x1=8, y1=7.5,
         lwd=lineWeights_Rates[5] * factor)    # 5  CowInner         -> CowOuterEast
  arrows(x0=4, y0=6.5, x1=2.5, y1=7.5,
         lwd=lineWeights_Rates[6] * factor)    # 6  CowInner         -> CowOuterWest
  
  arrows(x0=7, y0=7.5, x1=5.5, y1=6.5,
         lwd=lineWeights_Rates[5] * factor)    # 7  CowOuterEast     -> CowInner
  arrows(x0=7, y0=9, x1=3.5, y1=9,
         lwd=lineWeights_Rates[8] * factor)    # 8  CowOuterEast     -> CowOuterWest
  arrows(x0=9.5, y0=7.5, x1=9.5, y1=2.5,
         lwd=lineWeights_Rates[9] * factor)    # 9  CowOuterEast     -> BadgerOuterEast
  
  arrows(x0=3.5, y0=7.5, x1=5, y1=6.5,
         lwd=lineWeights_Rates[6] * factor)   # 10 CowOuterWest     -> CowInner
  arrows(x0=3.5, y0=8, x1=7, y1=8,
         lwd=lineWeights_Rates[8] * factor)   # 11 CowOuterWest     -> CowOuterEast
  arrows(x0=1, y0=7.5, x1=1, y1=2.5,
         lwd=lineWeights_Rates[9] * factor)   # 12 CowOuterWest     -> BadgerOuterWest
  
  arrows(x0=7, y0=2.5, x1=5.5, y1=3.5,
         lwd=lineWeights_Rates[3] * factor)   # 13 BadgerOuterEast  -> BadgerInner
  arrows(x0=8.5, y0=2.5, x1=8.5, y1=7.5,
         lwd=lineWeights_Rates[9] * factor)    # 14 BadgerOuterEast  -> CowOuterEast
  arrows(x0=7, y0=1, x1=3.5, y1=1,
         lwd=lineWeights_Rates[8] * factor)    # 15 BadgerOuterEast  -> BadgerOuterWest
  
  arrows(x0=3.5, y0=2.5, x1=5, y1=3.5,
         lwd=lineWeights_Rates[2] * factor)   # 16 BadgerOuterWest  -> BadgerInner
  arrows(x0=2, y0=2.5, x1=2, y1=7.5,
         lwd=lineWeights_Rates[9] * factor)  # 17 BadgerOuterWest  -> CowOuterWest
  arrows(x0=3.5, y0=2, x1=7, y1=2,
         lwd=lineWeights_Rates[8] * factor)   # 18 BadgerOuterWest  -> BadgerOuterEast

  text(x=c(4, 4, 0, 7.25, 7.25, 0), 
       y=c(4, 6, 8.5, 8.5, 1.5, 1.5), 
       labels=c("Inner Badgers", "Inner Cattle", "Outer West Cattle", "Outer East Cattle",
                "Outer East Badgers", "Outer West Badgers"),
       pos=4, col=c("red", "blue", "blue", "blue", rgb(1,0,0, 0.5), rgb(1,0,0, 0.5)), cex=cex)
  
}

plot8Deme <- function(path, popType, date, factor, cex){
  
  folderName <- paste("8Demes_", popType, "PopSizes_", date, "/", sep="")
  
  file <- paste(path, folderName,
                "8Demes_", popType, "PopSizes_", date, ".log", sep="")
  
  logTable <- read.table(file, header=TRUE, sep="\t")
  
  # Remove the burn-in
  burnIn <- round(0.1 * nrow(logTable), digits=0)
  logTable <- logTable[burnIn:nrow(logTable), ]
  
  #############################################
  # Calculate migration rates forward in time #
  #############################################
  
  logTable <- calculateForwardMigrationRates(logTable)
  
  ###################
  # Migration Rates #
  ###################
  
  # BadgerInnerEast   0     SAMPLED
  # BadgerInnerWest   1     SAMPLED
  # CowInnerEast      2     SAMPLED
  # CowInnerWest      3     SAMPLED
  # CowOuterEast      4     SAMPLED
  # CowOuterWest      5     SAMPLED
  # BadgerOuterEast   6     UNSAMPLED
  # BadgerOuterWest   7     UNSAMPLED
  
  # Set arrow direction
  code <- 2 # 2 for Backward in time (1 for forward)
  if(direction == "FORWARDS"){
    code <- 1
  }
  
  # Define line weights
  lineWeights_Rates <- c(
    mean(logTable[, "migModel.forwardRateMatrix_0_1"]), # 1  BadgerInnerEast  -> BadgerInnerWest
    mean(logTable[, "migModel.forwardRateMatrix_0_2"]), # 2  BadgerInnerEast  -> CowInnerEast
    mean(logTable[, "migModel.forwardRateMatrix_0_6"]), # 3  BadgerInnerEast  -> BadgerOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_1_0"]), # 4  BadgerInnerWest  -> BadgerInnerEast
    mean(logTable[, "migModel.forwardRateMatrix_1_3"]), # 5  BadgerInnerWest  -> CowInnerWest
    mean(logTable[, "migModel.forwardRateMatrix_1_7"]), # 6  BadgerInnerWest  -> BadgerOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_2_0"]), # 7  CowInnerEast     -> BadgerInnerEast
    mean(logTable[, "migModel.forwardRateMatrix_2_3"]), # 8  CowInnerEast     -> CowInnerWest
    mean(logTable[, "migModel.forwardRateMatrix_2_4"]), # 9  CowInnerEast     -> CowOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_3_1"]), # 10 CowInnerWest     -> BadgerInnerWest
    mean(logTable[, "migModel.forwardRateMatrix_3_2"]), # 11 CowInnerWest     -> CowInnerEast
    mean(logTable[, "migModel.forwardRateMatrix_3_5"]), # 12 CowInnerWest     -> CowOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_4_2"]), # 13 CowOuterEast     -> CowInnerEast
    mean(logTable[, "migModel.forwardRateMatrix_4_5"]), # 14 CowOuterEast     -> CowOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_4_6"]), # 15 CowOuterEast     -> BadgerOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_5_3"]), # 16 CowOuterWest     -> CowInnerWest
    mean(logTable[, "migModel.forwardRateMatrix_5_4"]), # 17 CowOuterWest     -> CowOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_5_7"]), # 18 CowOuterWest     -> BadgerOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_6_0"]), # 19 BadgerOuterEast  -> BadgerInnerEast
    mean(logTable[, "migModel.forwardRateMatrix_6_4"]), # 20 BadgerOuterEast  -> CowOuterEast
    mean(logTable[, "migModel.forwardRateMatrix_6_7"]), # 21 BadgerOuterEast  -> BadgerOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_7_1"]), # 22 BadgerOuterWest  -> BadgerInnerWest
    mean(logTable[, "migModel.forwardRateMatrix_7_5"]), # 23 BadgerOuterWest  -> CowOuterWest
    mean(logTable[, "migModel.forwardRateMatrix_7_6"])  # 24 BadgerOuterWest  -> BadgerOuterEast
  )
  
  ### Migration Rates
  plot(x=NULL, y=NULL, xlim=c(1,10), ylim=c(1,9.5),
       xaxt="n", yaxt="n", xlab="", ylab="", bty="n")
  
  arrows(x0=6, y0=3, x1=5, y1=3, 
         lwd=lineWeights_Rates[1] * factor)  # 1  BadgerInnerEast  -> BadgerInnerWest
  arrows(x0=7, y0=3, x1=7, y1=8, 
         lwd=lineWeights_Rates[2] * factor)  # 2  BadgerInnerEast  -> CowInnerEast
  arrows(x0=7.5, y0=3, x1=8, y1=4, 
         lwd=lineWeights_Rates[3] * factor)  # 3  BadgerInnerEast  -> BadgerOuterEast
 
  arrows(x0=5, y0=2, x1=6, y1=2, 
         lwd=lineWeights_Rates[4] * factor) # 4  BadgerInnerWest  -> BadgerInnerEast
  arrows(x0=4, y0=3, x1=4, y1=8, 
         lwd=lineWeights_Rates[5] * factor) # 5  BadgerInnerWest  -> CowInnerWest
  arrows(x0=3, y0=3, x1=2.5, y1=4, 
         lwd=lineWeights_Rates[6] * factor) # 6  BadgerInnerWest  -> BadgerOuterWest
   
  arrows(x0=6.5, y0=8, x1=6.5, y1=3, 
         lwd=lineWeights_Rates[7] * factor)  # 7  CowInnerEast     -> BadgerInnerEast
  arrows(x0=6, y0=9, x1=5, y1=9, 
         lwd=lineWeights_Rates[8] * factor)  # 8  CowInnerEast     -> CowInnerWest
  arrows(x0=7.5, y0=8, x1=8, y1=7, 
         lwd=lineWeights_Rates[9] * factor)  # 9  CowInnerEast     -> CowOuterEast
  
  arrows(x0=4.5, y0=8, x1=4.5, y1=3, 
         lwd=lineWeights_Rates[10] * factor)  # 10 CowInnerWest     -> BadgerInnerWest
  arrows(x0=5, y0=8, x1=6, y1=8, 
         lwd=lineWeights_Rates[11] * factor)  # 11 CowInnerWest     -> CowInnerEast
  arrows(x0=3, y0=8, x1=2.5, y1=7, 
         lwd=lineWeights_Rates[12] * factor)  # 12 CowInnerWest     -> CowOuterWest
  
  arrows(x0=9, y0=7, x1=8.5, y1=8, 
         lwd=lineWeights_Rates[13] * factor) # 13 CowOuterEast     -> CowInnerEast
  arrows(x0=7.5, y0=6, x1=3.5, y1=6, 
         lwd=lineWeights_Rates[14] * factor) # 14 CowOuterEast     -> CowOuterWest
  arrows(x0=9.5, y0=6, x1=9.5, y1=5, 
         lwd=lineWeights_Rates[15] * factor) # 15 CowOuterEast     -> BadgerOuterEast
  
  arrows(x0=1.5, y0=7, x1=2, y1=8, 
         lwd=lineWeights_Rates[16] * factor) # 16 CowOuterWest     -> CowInnerWest
  arrows(x0=3.5, y0=6.5, x1=7.5, y1=6.5, 
         lwd=lineWeights_Rates[17] * factor) # 17 CowOuterWest     -> CowOuterEast
  arrows(x0=1, y0=6, x1=1, y1=5, 
         lwd=lineWeights_Rates[18] * factor) # 18 CowOuterWest     -> BadgerOuterWest

  arrows(x0=9, y0=4, x1=8.5, y1=3, 
         lwd=lineWeights_Rates[19] * factor) # 19 BadgerOuterEast  -> BadgerInnerEast
  arrows(x0=8.5, y0=5, x1=8.5, y1=6, 
         lwd=lineWeights_Rates[20] * factor) # 20 BadgerOuterEast  -> CowOuterEast
  arrows(x0=7.5, y0=5, x1=3.5, y1=5, 
         lwd=lineWeights_Rates[21] * factor) # 21 BadgerOuterEast  -> BadgerOuterWest

  arrows(x0=1.5, y0=4, x1=2, y1=3, 
         lwd=lineWeights_Rates[22] * factor) # 22 BadgerOuterWest  -> BadgerInnerWest
  arrows(x0=2, y0=5, x1=2, y1=6, 
         lwd=lineWeights_Rates[23] * factor) # 23 BadgerOuterWest  -> CowOuterWest
  arrows(x0=3.5, y0=4.5, x1=7.5, y1=4.5, 
         lwd=lineWeights_Rates[24] * factor) # 24 BadgerOuterWest  -> BadgerOuterEast

  text(x=c(2, 6, 7.5, 7.5, 6, 2, 0.6, 0.60), 
       y=c(8.5, 8.5, 6.5, 4.5, 2.5, 2.5, 4.5, 6.5), 
       labels=c("Inner West Cattle", "Inner East Cattle", "Outer East Cattle",
                "Outer East Badgers", "Inner East Badgers", "Inner West Badgers", 
                "Outer West Badgers", "Outer West Cattle"),
       pos=4, col=c("blue", "blue", "blue", rgb(1,0,0, 0.5), "red", "red", rgb(1,0,0, 0.5), "blue"),
       cex=cex)
  
}

calculateForwardMigrationRates <- function(logTable){
  
  # Get the names of the backward in time migration rate estimates
  migrationRateCols <- colnames(logTable)[
    grepl(colnames(logTable), pattern = "migModel.rateMatrix")]
  
  # For each backward in time migration rate caculate the forward migration rate
  # FMR_ab = BMR_ba * (Nb / Na)
  #   MR: Migration rate (F - Forward, B - Backward)
  #   N: Effective population size
  #   Demes: a, b
  
  for(backwardMigrationRateCol in migrationRateCols){
    
    # Get the demes involved
    parts <- strsplit(backwardMigrationRateCol, split="_")[[1]]
    a <- parts[2]
    b <- parts[3]
    
    # Calculate forward rate
    forwardMigrationRateCol <- paste("migModel.forwardRateMatrix_", b, "_", a, sep="")
    logTable[, forwardMigrationRateCol] <-
      logTable[, backwardMigrationRateCol] *
      (logTable[, paste("migModel.popSize_", a, sep="")] * 
         logTable[, paste("migModel.popSize_", b, sep="")])
  }
  
  return(logTable)
}
