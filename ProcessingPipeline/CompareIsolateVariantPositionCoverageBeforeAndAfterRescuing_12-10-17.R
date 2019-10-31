############
# Set path #
############

path <- "/home/josephcrispell/Desktop/Research/EdgeArea_UK/vcfFiles/"

###############################################################
# Read in variant position coverage - before and after rescue #
###############################################################

# Read in the genome coverage file
file <- paste(path, "isolateCoverageSummary_DP-20_31-10-2019.txt", sep="")
coverage <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Read in the variant coverage file
file <- paste(path, "IsolateVariantPositionCoverage_FILTERED_31-10-2019.txt", sep="")
before <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

file <- paste(path, "IsolateVariantPositionCoverage_RESCUED_31-10-2019.txt", sep="")
after <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

#############################
# Plot the isolate coverage #
#############################

# Set the margins
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Open a pdf
file <- paste(substr(file, 1, nchar(file) - 4), ".pdf", sep="")
pdf(file)

table <- data.frame(Isolate=before$Isolate, stringsAsFactors=FALSE)
table$Before <- before$Coverage
table$After <- after$Coverage

values <- c(table$Before, table$After)

plot(table$Before, table$After, las=1, xlab="Before", ylab="After",
     main="Difference between before and after rescue",
     ylim=c(0,1), xlim=c(0,1), pch=19)

lines(x=range(values), y=range(values), lty=2, col="black")
for(row in 1:nrow(table)){
  
  lines(x=c(table[row, "Before"], table[row, "Before"]),
        y=c(table[row, "Before"], table[row, "After"]),
        lty=3, lwd=0.5, col="green")
}

plot(x=NULL, y=NULL, las=1, ylab="Coverage",
     main="Change in Coverage", xlab="", xaxt="n",
     ylim=c(0,1), xlim=c(0,1), bty="n")

for(row in 1:nrow(coverage)){
  
  lines(x=c(0.1, 0.5, 0.9),
        y=c(coverage[row, "PercentageCoverage"], 
            before[row, "Coverage"],
            after[row, "Coverage"]),
        type="o", lwd=0.5, col="black", pch=19)
}

axis(side=1, at=c(0.1, 0.5, 0.9), labels=c("Genome", "Before-VP", "After-VP"),
     tick="FALSE")

dev.off()