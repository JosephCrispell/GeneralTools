##### Preparation #####

# Set the path
path <- "/home/josephcrispell/Desktop/Research/Homoplasy/TimeTrial/"

#### Read in data ####

# Read in the table with the times
file <- paste(path, "timeTaken_27-04-18.csv", sep="")
table <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

#### Plot ####

# Open a PDF
file <- paste(path, "timeTaken_27-04-18.pdf", sep="")
pdf(file)

# Plot log of speeds
plot(x=table$N.sequences, y=log(table$R), type="o", las=1, bty="n",
     ylab="Time taken (seconds)", xlab="Number of sequences",
     main="Speed comparison between tools to identify homoplasies",
     yaxt="n", col="black", xaxt="n", ylim=c(-1, max(log(table$R))))
lines(x=table$N.sequences, y=log(table$java), type="o", col="red")
lines(x=table$N.sequences, y=log(table$treetime), type="o", col="blue")

# Add x axis
axis(side=1, at=seq(0, 500, 50))

# Add y axis (log scale)
at <- c(1, 10, 60, 150, 300, 600)
axis(side=2, at=c(-1,log(at)), labels=c(0, at), las=1)

# Add legend
legend("bottomright", legend=c("R", "Java", "TreeTime"), bty="n", text.col=c("black", "red", "blue"))

# Plot log of speeds
plot(x=table$N.sequences, y=table$R, type="o", las=1, bty="n",
     ylab="Time taken (seconds)", xlab="Number of sequences",
     main="Speed comparison between tools to identify homoplasies",
     col="black", ylim=c(0, max(table$R)), xaxt="n", yaxt="n")
lines(x=table$N.sequences, y=table$java, type="o", col="red")
lines(x=table$N.sequences, y=table$treetime, type="o", col="blue")

# Add x axis
axis(side=1, at=seq(50, 500, 50))

# Add y axis (log scale)
at <- c(0, 120, 240, 360, 480, 600)
axis(side=2, at=at, labels=at, las=1)

# Add legend
legend("topleft", legend=c("R", "Java", "TreeTime"), bty="n", text.col=c("black", "red", "blue"))

# Plot log of speeds
plot(x=table$N.sequences, y=table$java, type="o", las=1, bty="n",
     ylab="Time taken (seconds)", xlab="Number of sequences",
     main="Speed comparison between tools to identify homoplasies",
     col="red")
lines(x=table$N.sequences, y=table$treetime, type="o", col="blue")

# Add legend
legend("topleft", legend=c("Java", "TreeTime"), bty="n", text.col=c("red", "blue"))

dev.off()
      