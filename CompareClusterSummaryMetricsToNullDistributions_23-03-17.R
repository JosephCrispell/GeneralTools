# Create a path variable
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/InterSpeciesClusters/"

# Open the Cluster Summary table
file <- paste(path, "ClusterSummaryWithRandomNullDistributions_30-03-2017.txt", sep="")
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

file <- paste(path, "ClusterSummaryWithRandomNullDistributions_30-03-2017.pdf", sep="")

pdf(file, width=14, height=56)

par(mfrow=c(8, 2))
par(mar=c(5.1, 5.1, 1, 2.1)) # Bottom, Left, Top, Right

for(row in 1:nrow(summaryTable)){

  for(column in nameOrder){

    plotNullDistritionWithActual(summaryTable=summaryTable, nBreaks=15, 
                                 date=grepl(x=column, pattern="Date"),
                                 column=column, row=row, cexLab=2, cexAxis=1.5,
                                 cexLegend=2, cluster=row-1)
    
  }
}

dev.off()

#############
# FUNCTIONS #
#############

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