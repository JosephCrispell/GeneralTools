##### Preparation #####

# Set the path
path <- "/home/josephcrispell/storage/Research/Homoplasy/TimeTrial/"

#### Read in data ####

# Read in the table with the times
file <- paste(path, "timeTaken_05-12-19.csv", sep="")
table <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#### Plot ####

# Open a PDF
file <- paste(path, "timeTaken_05-12-19.pdf", sep="")
pdf(file)

# Create an empty plot on a log scale
plot(x=NULL, y=NULL, xlim=range(table$NSequences), ylim=range(table[, 4:ncol(table)]), col="white", las=1, bty="n",
     ylab="Time taken (seconds)", xlab="Number of sequences",
     main="Speed comparison between tools to identify homoplasies",
     xaxt="n")

# Add the mean values across 10 replicates and range
for(nSequences in seq(100, 1000, 50)){

  # Subset the table for the current value of sequences
  subset <- table[table$NSequences == nSequences, ]

  # R
  # points(x=nSequences, y=mean(subset$R), pch=19, col=rgb(0,0,0, 0.5))
  # points(x=c(nSequences, nSequences), y=range(subset$R), type="l", col=rgb(0,0,0, 0.5))

  # Java
  points(x=nSequences, y=mean(subset$Java), pch=19, col=rgb(1,0,0, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$Java), type="l", col=rgb(1,0,0, 0.5))

  # Phangorn
  # points(x=nSequences, y=mean(subset$Phangorn), pch=19, col=rgb(0,0,1, 0.5))
  # points(x=c(nSequences, nSequences), y=range(subset$Phangorn), type="l", col=rgb(0,0,1, 0.5))
  
  # TreeTime
  # points(x=nSequences, y=mean(subset$TreeTime), pch=19, col=rgb(0,1,1, 0.5))
  # points(x=c(nSequences, nSequences), y=range(subset$TreeTime), type="l", col=rgb(0,1,1, 0.5))
  
  # Java 2
  points(x=nSequences, y=mean(subset$Java2), pch=19, col=rgb(0,0,1, 0.5))
  points(x=c(nSequences, nSequences), y=range(subset$Java2), type="l", col=rgb(0,0,1, 0.5))
}

# Add x axis
axis(side=1, at=seq(100, 1000, 100))

# Add legend
# legend("topleft", legend=c("HomoplasyFinder (via R)", "HomoplasyFinder (via command line)", "Phangorn", "TreeTime"), 
#        bty="n", text.col=c("black", "red", "blue", "cyan"), cex=0.9)
legend("topleft", legend=c("OLD", "NEW"),
       bty="n", text.col=c("red", "blue"), cex=0.9)

dev.off()
