##### Preparation #####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Homoplasy/TimeTrial/"

#### Read in data ####

# Read in the table with the times
file <- paste(path, "timeTaken_17-08-18.csv", sep="")
table <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Plot ####

# Open a PDF
file <- paste(path, "timeTaken_17-08-18.pdf", sep="")
pdf(file)

# Create an empty plot on a log scale
plot(table$NSequences, y=table$TreeTime, col="white", las=1, bty="n",
     ylab="Time taken (seconds)", xlab="Number of sequences",
     main="Speed comparison between tools to identify homoplasies",
     xaxt="n", ylim=c(0, 50))

# Add the mean values across 10 replicates and range
for(nSequences in seq(100, 1000, 50)){

  # Subset the table for the current value of sequences
  subset <- table[table$NSequences == nSequences, ]

  # R
  points(x=nSequences, y=mean(subset$R), pch=19, col=rgb(0,0,0, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$R), type="l", col=rgb(0,0,0, 0.5))

  # Java
  points(x=nSequences, y=mean(subset$Java), pch=19, col=rgb(1,0,0, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$Java), type="l", col=rgb(1,0,0, 0.5))

  # Phangorn
  points(x=nSequences, y=mean(subset$Phangorn), pch=19, col=rgb(0,0,1, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$Phangorn), type="l", col=rgb(0,0,1, 0.5))
  
  # TreeTime
  points(x=nSequences, y=mean(subset$TreeTime), pch=19, col=rgb(0,1,1, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$TreeTime), type="l", col=rgb(0,1,1, 0.5))
}

# Add x axis
axis(side=1, at=seq(100, 1000, 100))

# Add legend
legend("topleft", legend=c("HomoplasyFinder (R)", "HomoplasyFinder (Java)", "Phangorn", "TreeTime"), 
       bty="n", text.col=c("black", "red", "blue", "cyan"))

dev.off()
