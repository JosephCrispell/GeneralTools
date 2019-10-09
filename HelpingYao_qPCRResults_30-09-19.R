# Open the qPCR data
fileName <- "/home/josephcrispell/Desktop/30_9_19_Sadie_linear_results.csv"
results <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")
results <- results[-1, ] # Remove the first row as it contains no data

# Log the data values
results[2:ncol(results)] <- log2(results[, 2:ncol(results)]) 

# Note the genes of interest
genes <- c("Mb2901", "Mb2902C", "SigA", "mpt70")

# Allow four plots in same window
par(mfrow=c(2,2))

# Create a plot for each gene
for(gene in genes){
  
  plotqPCRResults(results, geneOfInterest=gene, plotSE=FALSE)
}

#### FUNCTIONS ####

plotqPCRResults <- function(results, geneOfInterest, plotSE=TRUE){
  
  # Plot the data
  plot(x=NULL, y=NULL, xlim=c(0.5,3.5), ylim=range(results[, 2:ncol(results)], na.rm=TRUE),
       bty="n", las=1, xlab="", ylab="qPCR result (log2)", xaxt="n", main=geneOfInterest)
  
  # Add X axis
  axis(side=1, at=c(1,2,3), labels=c("wt", "ko", "c"))
  
  # Note point locations
  locations <- list("wt"=1, "ko"=2, "c"=3)
  pad <- c(-0.1, 0, 0.1)
  
  # Plot the data
  for(row in 1:nrow(results)){
    
    # Get the experiment type and replicate
    experiment <- strsplit(results[row, "Samples"], split="[0-9]", perl=TRUE)[[1]][1]
    replicate <- substr(results[row, "Samples"], nchar(results[row, "Samples"]), nchar(results[row, "Samples"]))
    
    # Note the X location
    xPosition <- locations[[experiment]] + pad[as.numeric(replicate)]
    
    # Plot the values
    mean <- results[row, paste0(geneOfInterest, ".CNRQ")]
    points(x=xPosition, y=mean, pch=19, col=rgb(0,0,0,0.5))
    
    # Plot the standard error if requested
    if(plotSE){
      standardError <- results[row, paste0(geneOfInterest, ".SE.CNRQ.")]
      points(x=c(xPosition, xPosition), y=c(mean-standardError, mean+standardError), type="l")
    }
  }
  
  # Run a wilcoxon rank test (mann-whitney) - complement versus wildtype
  valuesComplement <- results[grepl(results$Samples, pattern="c"), paste0(geneOfInterest, ".CNRQ")]
  valuesWildType <- results[grepl(results$Samples, pattern="wt"), paste0(geneOfInterest, ".CNRQ")]
  wilCoxRankTestComplement <- wilcox.test(x=valuesComplement, y=valuesWildType)
  text(x=locations[["c"]], y=mean(valuesComplement, na.rm=TRUE), labels=wilCoxRankTestComplement$p.value, pos=1)
  
  # Run a wilcoxon rank test (mann-whitney) - complement versus wildtype
  valuesKnockOut <- results[grepl(results$Samples, pattern="ko"), paste0(geneOfInterest, ".CNRQ")]
  wilCoxRankTestKnockout <- wilcox.test(x=valuesKnockOut, y=valuesWildType)
  text(x=locations[["ko"]], y=mean(valuesKnockOut, na.rm=TRUE), labels=wilCoxRankTestKnockout$p.value, pos=1)
}
