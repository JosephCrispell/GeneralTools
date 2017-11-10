##########
# Set up #
##########

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/BASTA/"

#####################
# Read in the files #
#####################

# Note deme structure to use
demeStructureDates <- list(
  "2Deme"="09-10-17",
  "3Deme-outerIsBoth"="09-10-17",
  "3Deme-outerIsCattle"="09-10-17",
  "4Deme"="09-10-17"#,
#  "6Deme-EastWest"="13-10-17",
#  "6Deme-NorthSouth"="13-10-17",
#  "8Deme-EastWest"="13-10-17",
#  "8Deme-NorthSouth"="13-10-17"
)

# Note the population size estimation options
popEstimationTypes <- c("varying", "equal")

# Note the clock model options
clockEstimateTypes <- "relaxed" # strict not used

# Store each of the log tables in a list
logTables <- readInBASTALogTables(demeStructureDates, popEstimationTypes, clockEstimateTypes,
                                  path)

####################
# Examine each run #
####################

# Get analysis names
analyses <- names(logTables)

# Examine each analysis
for(analysis in analyses){
  
  # Get the deme structure
  demeStructure <- 
  
  # Get the log table
  logTable <- logTables[[analysis]]
  
  # Plot the ESS values of each parameter estimated
  plotParameterEssValues(logTable, colNamesForESS)
  
  # Plot the posterior support for each deme as source
  plotPosteriorSupportForEachDemeAsRoot(logTable, demeStructure)
  
  # Produce a migration rate estimation figure - weight by rate flags
  # Diagrams designed with code = 2 (FORWARDS) in mind
  
}





arrowWeights <- list("0_1"=1, "1_0"=1, "0_2"=1, "2_0"=1, "1_2"=1,
                     "2_1"=1, "3_2"=1, "2_3"=1, "2_4"=1, "4_2"=1,
                     "3_5"=1, "5_3"=1, "4_5"=1, "5_4"=1, "0_4"=1,
                     "4_0"=1, "4_5"=1, "5_4"=1, "6_7"=1, "7_6"=1,
                     "4_6"=1, "6_4"=1, "5_7"=1, "7_5"=1, "0_6"=1,
                     "6_0"=1, "1_7"=1, "7_1"=1, "3_5"=1, "5_3"=1,
                     "4_6"=1, "6_4"=1)






#################################
# Plot Model comparison results #
#################################


######################################################################
# Produce single summary plot of rate estimation: BADGERS <-> CATTLE #
######################################################################

# Weight each model's input by model likelihood


#############
# FUNCTIONS #
#############

plot2Deme <- function(arrowWeights, code){

  demeStructure <- "2Deme"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.1, 0.9)
  y <- c(0.5, 0.5)
  text(x=x, y=y, 
       labels=demeNames,
       col=demeColours)
  
  # badger -> cow
  arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.15, y1=y[2]-0.15,
         code=code, lwd=arrowWeights[["1_0"]])
  # cow -> badger
  arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.15, y1=y[1]+0.15,
         code=code, lwd=arrowWeights[["0_1"]])
 
}

plot3Deme <- function(arrowWeights, code, demeStructure){
  
  demeStructure <- "3Deme-outerIsBoth"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Position the labels
  x <- c(0.1, 0.9, 0.5)
  y <- c(0.1, 0.1, 0.9)
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)
  
  # badger -> cow/cow-inner
  arrows(x0=x[1]+0.15, x1=x[2]-0.15, y0=y[1]-0.05, y1=y[2]-0.05,
         code=code, lwd=arrowWeights[["0_1"]]) 
  # cow/cow-inner -> badger
  arrows(x0=x[2]-0.15, x1=x[1]+0.15, y0=y[2]+0.05, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["1_0"]])
  
  # badger -> outer/cow-outer
  arrows(x0=x[1]-0.05, x1=x[3]-0.15, y0=y[1]+0.1, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["0_2"]])
  # outer/cow-outer -> badger
  arrows(x0=x[3]-0.05, x1=x[1]+0.05, y0=y[3]-0.1, y1=y[1]+0.1,
         code=code, lwd=arrowWeights[["2_0"]])
  
  # cow/cow-inner -> outer/cow-outer
  arrows(x0=x[2]+0.05, x1=x[3]+0.15, y0=y[2]+0.1, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["1_2"]])
  # outer/cow-outer -> cow/cow-inner
  arrows(x0=x[3]+0.05, x1=x[2]-0.05, y0=y[3]-0.1, y1=y[2]+0.1,
         code=code, lwd=arrowWeights[["2_1"]])
}

plot4Deme <- function(arrowWeights, code){
  
  demeStructure <- "4Deme"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Note the x and y positions
  x <- c(0.1, 0.1, 0.9, 0.9)
  y <- c(0.1, 0.9, 0.9, 0.1)
  
  # Add labels
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)

  # badger-inner -> cow-inner
  arrows(x0=x[1]-0.05, x1=x[2]-0.05, y0=y[1]+0.1, y1=y[2]-0.1,
         code=code, lwd=arrowWeights[["0_1"]]) 
  # cow-inner -> badger-inner
  arrows(x0=x[2]+0.05, x1=x[1]+0.05, y0=y[2]-0.1, y1=y[1]+0.1,
         code=code, lwd=arrowWeights[["0_1"]])
  
  # badger-outer -> cow-outer
  arrows(x0=x[4]+0.05, x1=x[3]+0.05, y0=y[4]+0.1, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["3_2"]]) 
  # cow-outer -> badger-outer
  arrows(x0=x[3]-0.05, x1=x[4]-0.05, y0=y[3]-0.1, y1=y[4]+0.1,
         code=code, lwd=arrowWeights[["2_3"]]) 
  
  # badger-inner -> badger-outer
  arrows(x0=x[1]+0.15, x1=x[4]-0.15, y0=y[1]-0.05, y1=y[4]-0.05,
         code=code, lwd=arrowWeights[["0_3"]]) 
  # badger-outer -> badger-inner
  arrows(x0=x[4]-0.15, x1=x[1]+0.15, y0=y[4]+0.05, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["3_0"]])
  
  # cow-inner -> cow-outer
  arrows(x0=x[2]+0.15, x1=x[3]-0.15, y0=y[2]+0.05, y1=y[3]+0.05,
         code=code, lwd=arrowWeights[["1_2"]]) 
  # cow-outer -> cow-inner
  arrows(x0=x[3]-0.15, x1=x[2]+0.15, y0=y[3]-0.05, y1=y[2]-0.05,
         code=code, lwd=arrowWeights[["2_1"]])
}

plot6DemeEW <- function(arrowWeights, code){
  
  demeStructure <- "6Deme-EastWest"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.5, 0.5, 0.9, 0.1, 0.9, 0.1)
  y <- c(0.35, 0.65, 0.9, 0.9, 0.1, 0.1)
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)
  
  # badger-inner -> cow-inner
  arrows(x0=x[1]-0.05, x1=x[2]-0.05, y0=y[1]+0.05, y1=y[2]-0.05,
         code=code, lwd=arrowWeights[["0_1"]]) 
  # cow-inner -> badger-inner
  arrows(x0=x[2]+0.05, x1=x[1]+0.05, y0=y[2]-0.05, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["1_0"]])
  
  # cow-outer-east -> badger-outer-east
  arrows(x0=x[3]-0.05, x1=x[5]-0.05, y0=y[3]-0.1, y1=y[5]+0.1,
         code=code, lwd=arrowWeights[["2_4"]]) 
  # badger-outer-east -> cow-outer-east
  arrows(x0=x[5]+0.05, x1=x[3]+0.05, y0=y[5]+0.1, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["4_2"]])
  
  # cow-outer-west -> badger-outer-west
  arrows(x0=x[4]+0.05, x1=x[6]+0.05, y0=y[4]-0.1, y1=y[6]+0.1,
         code=code, lwd=arrowWeights[["3_5"]]) 
  # badger-outer-west -> cow-outer-west
  arrows(x0=x[6]-0.05, x1=x[4]-0.05, y0=y[6]+0.1, y1=y[4]-0.1,
         code=code, lwd=arrowWeights[["5_3"]]) 
  
  # cow-outer-east -> cow-outer-west
  arrows(x0=x[3]-0.2, x1=x[4]+0.2, y0=y[3]-0.05, y1=y[4]-0.05,
         code=code, lwd=arrowWeights[["2_3"]]) 
  # cow-outer-west -> cow-outer-east 
  arrows(x0=x[4]+0.2, x1=x[3]-0.2, y0=y[4]+0.05, y1=y[3]+0.05,
         code=code, lwd=arrowWeights[["3_2"]])
  
  # badger-outer-east -> badger-outer-west
  arrows(x0=x[5]-0.2, x1=x[6]+0.2, y0=y[5]+0.05, y1=y[6]+0.05,
         code=code, lwd=arrowWeights[["4_5"]]) 
  # badger-outer-west -> badger-outer-east 
  arrows(x0=x[6]+0.2, x1=x[5]-0.2, y0=y[6]-0.05, y1=y[5]-0.05,
         code=code, lwd=arrowWeights[["5_4"]])
  
  # badger-inner -> badger-outer-east
  arrows(x0=x[1]+0.2, x1=x[5]-0.1, y0=y[1]-0.05, y1=y[5]+0.1,
         code=code, lwd=arrowWeights[["0_4"]])
  # badger-outer-east -> badger-inner
  arrows(x0=x[5]-0.2, x1=x[1]+0.1, y0=y[5]+0.1, y1=y[1]-0.05,
         code=code, lwd=arrowWeights[["0_4"]])
 
  # badger-inner -> badger-outer-west
  arrows(x0=x[1]-0.2, x1=x[6]+0.1, y0=y[1]-0.05, y1=y[6]+0.1,
         code=code, lwd=arrowWeights[["0_5"]])
  # badger-outer-west -> badger-inner
  arrows(x0=x[6]+0.2, x1=x[1]-0.1, y0=y[6]+0.1, y1=y[1]-0.05,
         code=code, lwd=arrowWeights[["5_0"]])
  
  # cow-inner -> cow-outer-east
  arrows(x0=x[2]+0.2, x1=x[3]-0.1, y0=y[2]+0.05, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["1_2"]])
  # cow-outer-east -> cow-inner
  arrows(x0=x[3]-0.2, x1=x[2]+0.1, y0=y[3]-0.1, y1=y[2]+0.05,
         code=code, lwd=arrowWeights[["2_1"]])

  # cow-inner -> cow-outer-west
  arrows(x0=x[2]-0.2, x1=x[4]+0.1, y0=y[2]+0.05, y1=y[4]-0.1,
         code=code, lwd=arrowWeights[["1_3"]])
  # cow-outer-west -> cow-inner
  arrows(x0=x[4]+0.2, x1=x[2]-0.1, y0=y[4]-0.1, y1=y[2]+0.05,
         code=code, lwd=arrowWeights[["3_1"]])
}

plot6DemeNS <- function(arrowWeights, code){
  
  demeStructure <- "6Deme-NorthSouth"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.35, 0.7, 0.9, 0.9, 0.1, 0.1)
  y <- c(0.5, 0.5, 0.9, 0.1, 0.9, 0.1)
  text(x=x, y=y, 
       labels=demeNames, col=demeColours)
  
  # badger-inner -> cow-inner
  arrows(x0=x[1]+0.1, x1=x[2]-0.1, y0=y[1]-0.05, y1=y[2]-0.05,
         code=code, lwd=arrowWeights[["0_1"]]) 
  # cow-inner -> badger-inner
  arrows(x0=x[2]-0.1, x1=x[1]+0.1, y0=y[2]+0.05, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["1_0"]])

  # cow-outer-north -> badger-outer-north
  arrows(x0=x[3]-0.2, x1=x[5]+0.2, y0=y[3]-0.05, y1=y[5]-0.05,
         code=code, lwd=arrowWeights[["2_4"]]) 
  # cow-outer-north -> badger-outer-north
  arrows(x0=x[5]+0.2, x1=x[3]-0.2, y0=y[5]+0.05, y1=y[3]+0.05,
         code=code, lwd=arrowWeights[["4_2"]])
  
  # cow-outer-south -> badger-outer-south
  arrows(x0=x[4]-0.2, x1=x[6]+0.2, y0=y[4]+0.05, y1=y[6]+0.05,
         code=code, lwd=arrowWeights[["3_5"]]) 
  # badger-outer-south -> cow-outer-south 
  arrows(x0=x[6]+0.2, x1=x[4]-0.2, y0=y[6]-0.05, y1=y[4]-0.05,
         code=code, lwd=arrowWeights[["5_3"]])
  
  # badger-outer-north -> badger-outer-south
  arrows(x0=x[5]+0.05, x1=x[6]+0.05, y0=y[5]-0.1, y1=y[6]+0.1,
         code=code, lwd=arrowWeights[["4_5"]]) 
  # badger-outer-south -> badger-outer-north
  arrows(x0=x[6]-0.05, x1=x[5]-0.05, y0=y[6]+0.1, y1=y[5]-0.1,
         code=code, lwd=arrowWeights[["5_4"]]) 
  
  # cow-outer-north -> cow-outer-south
  arrows(x0=x[3]-0.05, x1=x[4]-0.05, y0=y[3]-0.1, y1=y[4]+0.1,
         code=code, lwd=arrowWeights[["2_4"]]) 
  # cow-outer-south -> cow-outer-north
  arrows(x0=x[4]+0.05, x1=x[3]+0.05, y0=y[4]+0.1, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["3_2"]])
  
  # badger-inner -> badger-outer-north
  arrows(x0=x[1]-0.1, x1=x[5]+0.1, y0=y[1]+0.05, y1=y[5]-0.1,
         code=code, lwd=arrowWeights[["0_4"]])
  # badger-outer-north -> badger-inner
  arrows(x0=x[5]+0.2, x1=x[1], y0=y[5]-0.1, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["4_0"]])
  
  # badger-inner -> badger-outer-south
  arrows(x0=x[1]-0.1, x1=x[6]+0.1, y0=y[1]-0.05, y1=y[6]+0.1,
         code=code, lwd=arrowWeights[["0_5"]])
  # badger-outer-south -> badger-inner
  arrows(x0=x[6]+0.2, x1=x[1], y0=y[6]+0.1, y1=y[1]-0.05,
         code=code, lwd=arrowWeights[["5_0"]])
  
  # cow-inner -> cow-outer-north
  arrows(x0=x[2]+0.05, x1=x[3]-0.1, y0=y[2]+0.05, y1=y[3]-0.1,
         code=code, lwd=arrowWeights[["1_2"]])
  # cow-outer-north -> cow-inner
  arrows(x0=x[3]-0.2, x1=x[2]-0.05, y0=y[3]-0.1, y1=y[2]+0.05,
         code=code, lwd=arrowWeights[["2_1"]])
  
  # cow-inner -> cow-outer-south
  arrows(x0=x[2]+0.05, x1=x[4]-0.1, y0=y[2]-0.05, y1=y[4]+0.1,
         code=code, lwd=arrowWeights[["1_3"]])
  # cow-outer-south -> cow-inner
  arrows(x0=x[4]-0.2, x1=x[2]-0.05, y0=y[4]+0.1, y1=y[2]-0.05,
         code=code, lwd=arrowWeights[["3_1"]])

}

plot8DemeEW <- function(arrowWeights, code){
  
  demeStructure <- "8Deme-EastWest"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
  x <- c(0.75, 0.25, 0.75, 0.25, 0.9, 0.1, 0.9, 0.1)
  y <- c(0.1, 0.1, 0.9, 0.9, 0.65, 0.65, 0.35, 0.35)
  text(x=x, y=y, labels=demeNames, col=demeColours)
  
  # badger-inner-east -> cow-inner-east
  arrows(x0=x[1]-0.05, x1=x[3]-0.05, y0=y[1]+0.05, y1=y[3]-0.05,
         code=code, lwd=arrowWeights[["0_2"]]) 
  # cow-inner-east -> badger-inner-east
  arrows(x0=x[3]-0.125, x1=x[1]-0.125, y0=y[3]-0.05, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["2_0"]])
  
  # badger-inner-west -> cow-inner-west
  arrows(x0=x[2]+0.05, x1=x[4]+0.05, y0=y[2]+0.05, y1=y[4]-0.05,
         code=code, lwd=arrowWeights[["1_3"]]) 
  # cow-inner-west -> badger-inner-west
  arrows(x0=x[4]+0.125, x1=x[2]+0.125, y0=y[4]-0.05, y1=y[2]+0.05,
         code=code, lwd=arrowWeights[["3_1"]])
  
  # badger-inner-east -> badger-inner-west
  arrows(x0=x[1]-0.175, x1=x[2]+0.175, y0=y[1]+0.05, y1=y[2]+0.05,
         code=code, lwd=arrowWeights[["0_1"]]) 
  # badger-inner-west -> badger-inner-west
  arrows(x0=x[2]+0.175, x1=x[1]-0.175, y0=y[2]-0.05, y1=y[1]-0.05,
         code=code, lwd=arrowWeights[["1_0"]])
  
  # cow-inner-east -> cow-inner-west
  arrows(x0=x[3]-0.175, x1=x[4]+0.175, y0=y[3]-0.05, y1=y[4]-0.05,
         code=code, lwd=arrowWeights[["2_3"]]) 
  # cow-inner-west -> cow-inner-east
  arrows(x0=x[4]+0.175, x1=x[3]-0.175, y0=y[4]+0.05, y1=y[3]+0.05,
         code=code, lwd=arrowWeights[["3_2"]])
  
  # cow-outer-east -> cow-outer-west
  arrows(x0=x[5]-0.15, x1=x[6]+0.15, y0=y[5]-0.05, y1=y[6]-0.05,
         code=code, lwd=arrowWeights[["4_5"]]) 
  # cow-outer-west -> cow-outer-east
  arrows(x0=x[6]+0.15, x1=x[5]-0.15, y0=y[6]+0.05, y1=y[5]+0.05,
         code=code, lwd=arrowWeights[["5_4"]])
  
  # badger-outer-east -> badger-outer-west
  arrows(x0=x[7]-0.15, x1=x[8]+0.15, y0=y[7]+0.05, y1=y[8]+0.05,
         code=code, lwd=arrowWeights[["6_7"]]) 
  # badger-outer-west -> badger-outer-east
  arrows(x0=x[8]+0.15, x1=x[7]-0.15, y0=y[8]-0.05, y1=y[7]-0.05,
         code=code, lwd=arrowWeights[["7_6"]])
 
  # cow-outer-east -> badger-outer-east
  arrows(x0=x[5]-0.05, x1=x[7]-0.05, y0=y[5]-0.05, y1=y[7]+0.05,
         code=code, lwd=arrowWeights[["4_6"]]) 
  # badger-outer-east -> cow-outer-east
  arrows(x0=x[7]+0.05, x1=x[5]+0.05, y0=y[7]+0.05, y1=y[5]-0.05,
         code=code, lwd=arrowWeights[["6_4"]])
  
  # cow-outer-west -> badger-outer-west
  arrows(x0=x[6]+0.05, x1=x[8]+0.05, y0=y[6]-0.05, y1=y[8]+0.05,
         code=code, lwd=arrowWeights[["5_7"]]) 
  # badger-outer-west -> cow-outer-west
  arrows(x0=x[8]-0.05, x1=x[6]-0.05, y0=y[8]+0.05, y1=y[6]-0.05,
         code=code, lwd=arrowWeights[["7_5"]])
  
  # badger-inner-east -> badger-outer-east
  arrows(x0=x[1]+0.15, x1=x[7]+0.05, y0=y[1]+0.05, y1=y[7]-0.05,
         code=code, lwd=arrowWeights[["0_6"]]) 
  # badger-outer-east -> badger-inner-east
  arrows(x0=x[7]-0.05, x1=x[1]+0.05, y0=y[7]-0.05, y1=y[1]+0.05,
         code=code, lwd=arrowWeights[["6_0"]])
  
  # badger-inner-west -> badger-outer-west
  arrows(x0=x[2]-0.15, x1=x[8]-0.05, y0=y[2]+0.05, y1=y[8]-0.05,
         code=code, lwd=arrowWeights[["1_7"]]) 
  # badger-outer-west -> badger-inner-west
  arrows(x0=x[8]+0.05, x1=x[2]-0.05, y0=y[8]-0.05, y1=y[2]+0.05,
         code=code, lwd=arrowWeights[["7_1"]])
  
  # cow-inner-east -> cow-outer-east
  arrows(x0=x[3]+0.15, x1=x[5]+0.05, y0=y[3]-0.05, y1=y[5]+0.05,
         code=code, lwd=arrowWeights[["2_4"]]) 
  # cow-outer-east -> cow-inner-east
  arrows(x0=x[5]-0.05, x1=x[3]+0.05, y0=y[5]+0.05, y1=y[3]-0.05,
         code=code, lwd=arrowWeights[["4_2"]])
  
  # cow-inner-west -> cow-outer-west
  arrows(x0=x[4]-0.15, x1=x[6]-0.05, y0=y[4]-0.05, y1=y[6]+0.05,
         code=code, lwd=arrowWeights[["3_5"]]) 
  # cow-outer-west -> cow-inner-west
  arrows(x0=x[6]+0.05, x1=x[4]-0.05, y0=y[6]+0.05, y1=y[4]-0.05,
         code=code, lwd=arrowWeights[["5_3"]])
  
}

plot8DemeNS <- function(arrowWeights, code){
  
  demeStructure <- "8Deme-NorthSouth"
  
  # Create empty plot
  createEmptyPlot()
  
  # Get deme names and assign colours
  demeNames <- getDemeNamesForDemeStructure(demeStructure)
  demeColours <- rep("black", length(demeNames))
  demeColours[grepl(demeNames, pattern="badger")] <- "red"
  demeColours[grepl(demeNames, pattern="cow")] <- "blue"
  
  # Add labels
#  x <- c(0.25, 0.25, 0.75, 0.75, 0.9,  0.9,  0.1,  0.1)
#  y <- c(0.9,  0.1,  0.9,  0.1,  0.65, 0.35, 0.65, 0.35)

#  x <- c(0.1, 0.1, 0.9, 0.9, 0.75,  0.75,  0.25,  0.25)
#  y <- c(0.65,  0.35,  0.65,  0.35,  0.9, 0.1, 0.9, 0.1)
  
  x <- c(0.25, 0.25, 0.75, 0.75, 0.9,  0.9,  0.1,  0.1)
  y <- c(0.65,  0.35,  0.65,  0.35,  0.9, 0.1, 0.9, 0.1)
  
  text(x=x, y=y, labels=demeNames, col=demeColours)
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

createEmptyPlot <- function(){
  par(mar=c(0,0,0,0))
  plot(x=NULL, y=NULL, xlim=c(0,1), ylim=c(0,1), xaxt="n", yaxt="n", bty="n",
       ylab="", xlab="")
}

plotPosteriorSupportForEachDemeAsRoot <- function(logTable, demeStructure){

  barplot(table(logTable$treePrior.rootColor), yaxt="n", names=getDemeNamesForDemeStructure(demeStructure),
          main="Assignment of Demes to Root State")
}

plotParameterESSValues <- function(logTable, colNamesToPlot){
  
  essValues <- c()

  for(i in 1:length(colsToPlot)){
    
    essValues[i] <- calculateEffectiveSampleSize(logTable[, colNamestoPlot[i]])
  }
  
  par(mar=c(0,11,2,0.5)) # bottom, left, top, right
  barplot(essValues, las=1, names=colNamesToPlot, horiz=TRUE, xaxt="n", cex.names=0.9,
          main="Parameter Effective Sample Sizes")
  abline(v=100, lty=2, col="red")
  abline(v=1000, lty=2, col="blue")
  text(x=c(250, 1200), y=c(0, 0), labels=c("100", "1000"), cex=0.5, col=c("red", "blue"))
  
  # Reset Margins
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

calculateEffectiveSampleSize <- function(posteriorSample){
  
  # Calculate the mean of the sample
  sampleMean <- mean(posteriorSample)
  nSamples <- length(posteriorSample)
  
  # Create an array to store the correlation values between a[index] vs. a[index + lag]
  autoCorrelationValues <- c()
  
  for(lag in 1:(nSamples - 1)){
    
    # Calculate the auto-correlation value for the current lag period
    # Taken from: http://www.itl.nist.gov/div898/handbook/eda/section3/autocopl.htm
    a <- posteriorSample[1:(nSamples - lag)] - sampleMean
    b <- posteriorSample[(1:(nSamples - lag) + lag)] - sampleMean
    
    Ch <- (1/nSamples) * sum(a * b)
    C0 <- sum((posteriorSample - sampleMean)^2) / nSamples
    autoCorrelationValues[lag] <- Ch / C0
    
    # Check if auto-correlation has dropped below zero
    if(lag != 1 && 
       (autoCorrelationValues[lag - 1] > 0 && autoCorrelationValues[lag] <= 0 || 
        autoCorrelationValues[lag - 1] < 0 && autoCorrelationValues[lag] >= 0)){
      break;
    }
    
    # Monitor progress
    #if(lag %% 1000 == 0){
    #  cat(paste("Calculate correlation based upon gap of ", lag, ". Max gap = ", (nSamples - 1), "\n", sep=""))
    #}
  }
  
  # Calculate the Effective Sample Size
  # Taken from: http://people.duke.edu/~ccc14/sta-663-2016/16C_PyMC3.html
  ess <- nSamples / (1 + 2 * (sum(autoCorrelationValues)))
  
  return(ess);
}

getDemeNamesForDemeStructure <- function(demeStructure){

  # Species-InOrOut-Location
  # badger  inner   east
  # cow     outer   west
  #                 north
  #                 south
  demeNames <- list(
    "2Deme"=c("badger", "cow"),
    
    "3Deme-outerIsBoth"=c("badger-inner", "cow-inner", "outer"),
    
    "3Deme-outerIsCattle"=c("badger-inner", "cow-inner", "cow"),
    
    "4Deme"=c("badger-inner", "cow-inner", "cow-outer", "badger-outer"),
    
    "6Deme-EastWest"=c("badger-inner", "cow-inner", "cow-outer-east", "cow-outer-west", "badger-outer-east", "badger-outer-west"),
    
    "6Deme-NorthSouth"=c("badger-inner", "cow-inner", "cow-outer-north", "cow-outer-south", "badger-outer-north", "badger-outer-south"),
    
    "8Deme-EastWest"=c("badger-inner-east", "badger-inner-west", "cow-inner-east", "cow-inner-west",
                       "cow-outer-east", "cow-outer-west", "badger-outer-east", "badger-outer-west"),
    
    "8Deme-NorthSouth"=c("badger-inner-north", "badger-inner-south", "cow-inner-north", "cow-inner-south",
                         "cow-outer-north", "cow-outer-south", "badger-outer-north", "badger-outer-south")
  )
  
  return(demeNames[[demeStructure]])
}

readInBASTALogTables <- function(demeStructureDates, popEstimationTypes, clockEstimateTypes,
                                 path){
  
  # Store each of the log tables in a list
  logTables <- list()
  
  for(demeStructure in names(demeStructureDates)){
    
    for(popEstimationType in popEstimationTypes){
      
      for(clockEstimationType in clockEstimateTypes){
        
        # Build run defining prefix
        prefix <- paste(demeStructure, "_", popEstimationType, "_",  clockEstimationType, "_",
                        demeStructureDates[[demeStructure]], sep="")
        
        # Print progress information
        cat(paste("\rReading: ", prefix, ".log                    ", sep=""))
        
        # Create file name
        file <- paste(path, prefix, "/", prefix, ".log", sep="")
        
        # Read in the file as table
        logTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        
        # Replace "N" in table with NAs
        logTable[logTable == "N"] <- NA
        
        # Remove the burn-in
        burnIn <- round(0.1 * nrow(logTable), digits=0)
        logTable <- logTable[burnIn:nrow(logTable), ]
        
        # Calculate the forward rates
        logTable <- calculateForwardMigrationRates(logTable)
        
        # Store the table
        logTables[[prefix]] <- logTable
      }
    }
  }
  cat("\rFinished reading in log tables...                              ")
  
  return(logTables)
}
