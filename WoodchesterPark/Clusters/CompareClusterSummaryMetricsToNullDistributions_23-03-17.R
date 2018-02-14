# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/InterSpeciesClusters/"

# Open the Cluster Summary table
file <- paste(path, "ClusterSummaryWithRandomNullDistributions_13-02-2018.txt", sep="")
summaryTable <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Define the order of the plots
nameOrder <- c(
  "NUnSampledDetectedBadgers", "NUnSampledDetectedCattle",
  "NNegativeBadgers", "NNegativeCattle",
  "EarliestDetectionDateSampledBadgers", "EarliestDetectionDateSampledCattle",
  "EarliestDetectionDateUnSampledBadgers", "EarliestDetectionDateUnSampledCattle",
  "MeanShortestPathLengthBetweenSampledGroups", "ProportionShortestPathsBetweenSampledGroupsThatExist",
  "MeanShortestPathLengthBetweenSampledHerds", "ProportionShortestPathsBetweenSampledHerdsThatExist",
  "NumberSampledGroups","NumberSampledHerds",
  
  "NUnSampledInconclusive",
  "MeanSpatialDist"
)

# file <- paste(path, "ClusterSummaryWithRandomNullDistributions_30-03-2017.pdf", sep="")
# 
# pdf(file, width=14, height=56)
# 
# par(mfrow=c(8, 2))
# par(mar=c(5.1, 5.1, 1, 2.1)) # Bottom, Left, Top, Right
# 
# for(row in 1:nrow(summaryTable)){
# 
#   for(column in nameOrder){
# 
#     plotNullDistritionWithActual(summaryTable=summaryTable, nBreaks=15, 
#                                  date=grepl(x=column, pattern="Date"),
#                                  column=column, row=row, cexLab=2, cexAxis=1.5,
#                                  cexLegend=2, cluster=row-1)
#     
#   }
# }
# 
# dev.off()

#######################################################
# Build an output table with Null Distribution bounds #
#######################################################

# Create the output table
fullNames <- list(
  "MeanDistToRefBadgers" = "Mean number of SNPs from clade root for badger sequences",
  "MeanDistToRefCattle" = "Mean number of SNPs from clade root for cattle sequences",
  "MeanSeqQualBadgers" = "Mean genome coverage of badger sequences",
  "MeanSeqQualCattle" = "Mean genome coverage of cattle sequences",
  "NBadgersSampled" = "Number of badgers sampled",
  "NCattleSampled" = "Number of cattle sampled",
  "NUnSampledDetectedBadgers" = "Number of in-contact badgers that tested positive",
  "NUnSampledDetectedCattle" = "Number of in-contact cattle that tested positive",
  "NUnSampledInconclusive" = "Number of in-contact cattle whose test repsonse was inconclusive",
  "NNegativeBadgers" = "Number of in-contact badgers that NEVER tested positive",
  "NNegativeCattle" = "Number of in-contact cattle that NEVER tested positive",
  "EarliestDetectionDateSampledBadgers" = "Earliest date that a sampled badger tested positive",
  "EarliestDetectionDateSampledCattle" = "Earliest date that a sampled cow tested positive",
  "EarliestDetectionDateUnSampledBadgers" = "Earliest date that an in-contact badger tested positive",
  "EarliestDetectionDateUnSampledCattle" = "Earliest date that an in-contact cow tested positive",
  "MeanSpatialDist" = "Mean spatial distance (KM) from the sampled herds to Woodchester Park",
  "MeanShortestPathLengthBetweenSampledGroups" = "Mean length of shortest paths between sampled social groups",
  "MeanShortestPathLengthBetweenSampledHerds" = "Mean length of shortest paths between sampled herds",
  "ProportionShortestPathsBetweenSampledGroupsThatExist" = "Proportion of shortest paths that exist between sampled social groups",
  "ProportionShortestPathsBetweenSampledHerdsThatExist" = "Proportion of shortest paths that exist between herds",
  "NumberSampledGroups" = "Number of social groups sampled",
  "NumberSampledHerds" = "Number of herds sampled"
)

colNames <- c(
  "Cluster-0",
  "Cluster-1",
  "Cluster-2",
  "Cluster-3"
)

outputTable <- buildOutputTable(summaryTable, fullNames, colNames)

# Replace the row names
rownames(outputTable) <- replaceFromList(rownames(outputTable), fullNames)

# Save the output table as csv
file <- paste(path, "ClusterSummaryWithBoundsOfNulls_13-02-2018.txt", sep="")
write.table(file, x=outputTable, quote=FALSE, row.names=TRUE, sep="\t")


#############
# FUNCTIONS #
#############

replaceFromList <- function(array, list){
  
  output <- c()
  for(i in 1:length(array)){
    output[i] <- list[[array[i]]]
  }
  
  return(output)
}

buildOutputTable <- function(summaryTable, fullNames, colNames){

  outputTable <- data.frame(matrix(nrow=length(names(fullNames)), ncol=length(colNames)))
  colnames(outputTable) <- colNames
  rownames(outputTable) <- names(fullNames)
  
  # Fill in the output table
  for(cluster in 0:(nrow(summaryTable)-1)){
    for(metric in names(fullNames)){
      
      # Check if current cell is dealing with a date
      if(grepl(x=metric, pattern="Date") == TRUE){
        
        # Split the current cell into its parts
        values <- as.Date(strsplit(summaryTable[cluster + 1, metric], split=",")[[1]],
                          format="%d-%m-%Y")
        
        # Find the bounds of the null distribution
        bounds <- quantile(values[-1], probs=c(0.025, 0.975), type=1, na.rm=TRUE) # type 1 - algorithm to handle discrete distribution
        
        # Put cell info into output table
        outputTable[metric, cluster + 1] <- paste(values[1], " (", bounds[[1]], ", ", bounds[[2]], ")", sep="")
        
        # Check if cell has no null distribution
      }else if(grepl(x=metric, pattern="DistToRef") == TRUE || 
               grepl(x=metric, pattern="SeqQual") == TRUE){
        
        # Put cell info into output table
        outputTable[metric, cluster + 1] <- round(summaryTable[cluster + 1, metric], digits=2)
        
        # Don't examine null distribution for N sampled metrics - all the same
      }else if(metric == "NBadgersSampled" || metric == "NCattleSampled"){
        
        # Split the current cell into its parts
        values <- as.numeric(strsplit(summaryTable[cluster + 1, metric], split=",")[[1]])
        
        # Put cell info into output table
        outputTable[metric, cluster + 1] <- values[1]
        
        # Summarise null distribution of numeric data
      }else{
        
        # Split the current cell into its parts
        values <- as.numeric(strsplit(summaryTable[cluster + 1, metric], split=",")[[1]])
        
        # If spatial dists - divide by 1000 - convert to km
        if(metric == "MeanSpatialDist"){
          values <- values / 1000
        }
        
        # Find the bounds of the null distribution
        if(length(values) > 1){
          bounds <- round(quantile(values[-1], probs=c(0.025, 0.975)), digits=2)
        }
        
        # Put cell info into output table
        outputTable[metric, cluster + 1] <- paste(round(values[1], digits=2),
                                                  " (", bounds[[1]], ", ", bounds[[2]], ")", sep="")
      }
    }
  }
  
  return(outputTable)
}

plotNullDistritionWithActual <- function(summaryTable, column, row, nBreaks, date,
                                         cexLab, cexAxis, cexLegend, cluster){
  
  # Split the cell into its values
  if(date == FALSE){
    values <- as.numeric(strsplit(summaryTable[row, column], split=",")[[1]])
    
    h <- hist(values[-1], breaks=nBreaks, las=1, freq=TRUE, 
              xlab=column, main="", cex.lab=cexLab, cex.axis=cexAxis,
              yaxt="n")
  }else{
    values <- as.Date(strsplit(summaryTable[row, column], split=",")[[1]],
                      format="%d-%m-%Y")
    
    h <- hist(values[-1], breaks="years", las=1, freq=TRUE, 
              xlab=column, main="", cex.lab=cexLab, cex.axis=cexAxis,
              format="%d-%m-%Y", xlim=as.Date(c("1985-01-01", "2015-01-01")),
              xaxt="n")
    axis.Date(side=1, at=seq(from=as.Date("1985-01-01"), to=as.Date("2015-01-01"), by="year"),
         format="%d-%m-%Y", cex.axis=cexAxis)
  }
  
  # Add a line to show the actual value
  if(cluster != 5){
    lines(c(values[1],values[1]), c(0, max(h$counts)), col="blue", lwd=3, xpd=TRUE)
  }
  
  
  legend("topright", legend=paste("Cluster:", cluster), bty="n", cex=cexLegend)
  
}