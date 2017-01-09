
# Create two random samples from a normal distribution
sampleA <- rnorm(50, mean=5) # Draw 50 random samples from a normal distribution with a mean of 5
sampleB <- rnorm(50, mean=10) # Draw 50 random samples from a normal distribution with a mean of 5

## Plot them in a histogram

# Define the X-axis limits
xLim <- c(min(sampleA, sampleB), max(sampleA, sampleB))

# Define the bins
# Define a sequence of numbers that are the bounds of bins for histogram
breaks <- seq(from=xLim[0], to=xLim[1] + 0.1, by=0.25) 

# Create the histograms but don't plot them yet
aHist <- hist(sampleA, breaks=breaks, plot=FALSE)
bHist <- hist(sampleB, breaks=breaks, plot=FALSE)

# Define the Y-axis limits
yLim <- c(0, max(aHist$counts, bHist$counts))

# Define the colours
aColour <- rgb(1,0,0, 0.5) # Red, Green, Blue, Transparency     RED
bColour <- rgb(0,0,1, 0.5) # Red, Green, Blue, Transparency     BLUE

# Plot the histograms
plot(aHist, col=aColour)
plot(bHist, col=bColour, add=TRUE)
