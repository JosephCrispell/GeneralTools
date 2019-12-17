############
# Set path #
############

path <- "/home/josephcrispell/storage/Research/RepublicOfIreland/Mbovis/Monaghan/Fastqs_16-12-19/"

###################################
# Read in isolate mapping summary #
###################################

# Note the file name
file <- paste(path, "isolateMappingSummary_16-12-19.txt", sep="")

# Read in the table
mappingSummary <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

# Parse the Isolate column
mappingSummary$Isolate <- parseIsolateColumn(mappingSummary$Isolate)

################################
# Plot the mapping information #
################################

# Calculate mapping proportion
mappingSummary$ProportionMapped <- mappingSummary$NumberMappedReads / 
  (mappingSummary$NumberMappedReads + mappingSummary$NumberUnmappedReads)

# Open a pdf
file <- paste(substr(file, 1, nchar(file) - 4), ".pdf", sep="")
pdf(file)

plot(y=mappingSummary$ProportionMapped,
     x=mappingSummary$NumberMappedReads + mappingSummary$NumberUnmappedReads,
     las=1, ylab="Proportion", main="Proportion Reads Mapped to M. bovis genome",
     xlab="Total Number of Reads", pch=16, col=rgb(0,0,0, 0.5), bty="n")

text(y=mappingSummary$NumberMappedReads / 
       (mappingSummary$NumberMappedReads + mappingSummary$NumberUnmappedReads),
     x=mappingSummary$NumberMappedReads + mappingSummary$NumberUnmappedReads,
     labels = mappingSummary$Isolate, cex=0.5, xpd=TRUE,
     col=ifelse(mappingSummary$ProportionMapped < 0.9, rgb(1,0,0, 1), rgb(0,0,0, 0)))

dev.off()

#############
# FUNCTIONS #
#############

parseIsolateColumn <- function(column){
  
  ids <- c()
  for(i in 1:length(column)){
    
    ids[i] <- strsplit(column[i], split="_")[[1]][1]
  }
  
  return(ids)
}
