##### Preparation #####

# Set the path
#path <- "/home/josephcrispell/Desktop/Research/Homoplasy/TimeTrial/"
#path <- "/home/josephcrispell/Desktop/Research/Homoplasy/TimeTrial/"
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Homoplasy/TimeTrial/"

#### Read in data ####

# Read in the table with the times
file <- paste(path, "timeTaken_11-05-18.csv", sep="")
table <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#### Plot ####

# Open a PDF
file <- paste(path, "timeTaken_11-05-18.pdf", sep="")
pdf(file)

# Create an empty plot on a log scale
plot(table$NSequences, y=table$R, col="white", las=1, bty="n",
     ylab="Time taken (seconds)", xlab="Number of sequences",
     main="Speed comparison between tools to identify homoplasies",
     xaxt="n", ylim=c(-1.1, max(table$R)))

# Add the mean values across 10 replicates and range
for(nSequences in seq(50, 500, 50)){

  # Subset the table for the current value of sequences
  subset <- table[table$NSequences == nSequences, ]

  # R
  points(x=nSequences, y=mean(subset$R), pch=19, col=rgb(0,0,0, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$R), type="l", col=rgb(0,0,0, 0.5))

  # Java
  points(x=nSequences, y=mean(subset$Java), pch=19, col=rgb(1,0,0, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$Java), type="l", col=rgb(1,0,0, 0.5))

  # TreeTime
  points(x=nSequences, y=mean(subset$TreeTime), pch=19, col=rgb(0,0,1, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$TreeTime), type="l", col=rgb(0,0,1, 0.5))
}

# Add x axis
axis(side=1, at=seq(50, 500, 50))

# Add legend
legend("bottomright", legend=c("R", "Java", "TreeTime"), bty="n", text.col=c("black", "red", "blue"))

dev.off()
