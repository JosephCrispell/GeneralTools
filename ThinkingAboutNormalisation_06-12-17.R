######################################
# A script to describe normalisation #
######################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/"

# Open a PNG
file <- paste(path, "ThinkingaboutNormalisation_06-12-17.png", sep="")
png(file, width=960, height=960)

# Some values 
values <- c(2,3,4,3,5,6,8,9)

# Add 10 to all values
values <- values + 10

# Split the plotting window
par(mfrow=c(2,4))

# Plot those values
plot(x=rep(1, length(values)), y=values,
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), main="True values",
     cex.main=2, cex.axis=2)

# Divide values by max
dividedByMax <- values/max(values)

# Plot those values
plot(x=rep(1, length(dividedByMax)), y=dividedByMax,
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), main="Divided by max(values)",
     cex.main=2, cex.axis=2)

# Normalise the values
normalised <- (values - min(values)) / (max(values) - min(values))

# Plot the normalised values
plot(x=rep(1, length(normalised)), y=normalised,
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), 
     main="Normalised",
     cex.main=2, cex.axis=2)

# Re-set scale of values
minToSet <- 5
maxToSet <- 50
rescaled <- setRangeOfValues(values, minToSet, maxToSet)

# Plot the rescaled values
plot(x=rep(1, length(rescaled)), y=rescaled,
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), 
     main=paste("Re-scaled between", minToSet, "and", maxToSet),
     cex.main=2, cex.axis=2)

# Plot those values
plot(x=rep(1, length(values)), y=values,
     ylim=c(0, max(values)),
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), main="True values",
     cex.main=2, cex.axis=2)

# Divide values by max
dividedByMax <- values/max(values)

# Plot those values
plot(x=rep(1, length(dividedByMax)), y=dividedByMax,
     ylim=c(0, max(dividedByMax)),
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), main="Divided by max(values)",
     cex.main=2, cex.axis=2)

# Normalise the values
normalised <- (values - min(values)) / (max(values) - min(values))

# Plot the normalised values
plot(x=rep(1, length(normalised)), y=normalised,
     ylim=c(0, max(normalised)),
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), 
     main="Normalised",
     cex.main=2, cex.axis=2)

# Re-set scale of values
minToSet <- 5
maxToSet <- 50
rescaled <- setRangeOfValues(values, minToSet, maxToSet)

# Plot the rescaled values
plot(x=rep(1, length(rescaled)), y=rescaled,
     ylim=c(0, max(rescaled)),
     xaxt="n", ylab="", xlab="", las=1, bty="n",
     pch=19, cex=4, col=rgb(0,0,0, 0.75), 
     main=paste("Re-scaled between", minToSet, "and", maxToSet),
     cex.main=2, cex.axis=2)

dev.off()

#######################################################################################
# Normalisation                                                                       #
# - Only provides difference between a and b relative to distances of a to all others #
# - Destroys vector with only two elements - where only interesting information is    #
#   the difference between a and b                                                    #
#######################################################################################


#############
# FUNCTIONS #
#############

setRangeOfValues <- function(values, minToSet=0, maxToSet=1){
  
  # Took the normalisation equation from:
  #   https://en.wikipedia.org/wiki/Feature_scaling
  
  # Create a vector to store the normalised values
  output <- c()
  
  # Calculate the min and max of the input values
  range <- range(values)
  
  # Normalise each of the input values and set on scale of the input min and max scales
  for(i in 1:length(values)){
    
    # Normalise between 0 and maxToSet - minToSet
    output[i] <- ((values[i] - range[1]) / (range[2] - range[1])) * (maxToSet - minToSet)
    
    # Push distribution up by minToSet
    output[i] <- output[i] + minToSet
  }
  
  return(output)
}
