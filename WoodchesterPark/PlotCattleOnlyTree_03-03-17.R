#############
# Libraries #
#############

library(ape)
library(geiger)
library(randomForest)

########
# Path #
########

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

#################################
# Open the ML Phylogenetic Tree #
#################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Convert branch lengths into SNPs
tree$edge.length <- tree$edge.length * 9464

# Remove all the badger isolates
tree <- removeBadgersFromTree(tree)

# Index the tips
indexedTips <- indexArray(tree$tip.label)

###################
# Get the Eartags #
###################

# Get the sampling information
file <- paste(path, "IsolateData/",
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
samplingInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

# Get the eartags for each tip label
tipInfo <- getTipInfo(indexedTips, samplingInfo, c("Rawtag", "Mapx", "Mapy"))

##########################
# Plot the sampled herds #
##########################

# Open a PDF
file <- paste(path, "CattleOnlyTreeWithEartags_03-03-17.pdf", sep="")
pdf(file, height=14, width=14)

# Set the margins
par(mar=c(5.1,5.1,4.1,2.1))

# Calculate the ranges of the x and y axes
minX <- min(tipInfo$Mapx, na.rm=TRUE)
minY <- min(tipInfo$Mapy, na.rm=TRUE)
axisSize <- 25000

plot(x=tipInfo$Mapx, y=tipInfo$Mapy, pch=19, col=rgb(0,0,0, 0.25),
     xlim=c(minX, minX + axisSize), ylim=c(minY, minY + axisSize),
     xaxt="n", yaxt="n", asp=1, bty="n", 
     xlab=paste(axisSize / 1000, " km", sep=""),
     ylab=paste(axisSize / 1000, " km", sep=""),
     main="Relative locations of Sampled Herds", cex.main=3, cex=4, cex.lab=2.5)

#################
# Plot the Tree #
#################

# Set the margins
par(mar=c(0,0,0,0))

# Change tip labels to eartags
treeEartags <- tree
treeEartags$tip.label <- tipInfo$Rawtag
treeEartags$tip.label[length(treeEartags$tip.label)] <- "Reference"

# Plot the phylogenetic tree with defined clades
plot.phylo(treeEartags, show.tip.label=TRUE, "fan",
           edge.color="black", edge.width=3, label.offset=1)

# Add a Scale bar
points(x=c(95, 145), y=c(-100, -100), type="l", lwd=3, col="grey")
text(x=120, y=-105, labels="50 SNPs", cex=1)

plot.new()
text(x=rep(0.5, length(tipInfo$Rawtag)), y=seq(1, 0, -1/length(tipInfo$Rawtag)),
     labels=tipInfo$Rawtag)

dev.off()



#############
# FUNCTIONS #
#############

getTipInfo <- function(indexedTips, samplingInfo, columns){
  
  tipInfo <- data.frame(TipLabel=names(indexedTips), stringsAsFactors=FALSE)
  
  # Add the other columns
  tipInfo[, columns] <- rep(NA, nrow(tipInfo))
  
  # Examine each row of the sampling info table
  for(row in 1:nrow(samplingInfo)){
    
    # Is this information associated with an isolate we're interested in?
    if(is.null(indexedTips[[samplingInfo[row, "StrainId"]]]) == FALSE){
      
      # Store its eartag
      for(col in columns){
        tipInfo[indexedTips[[samplingInfo[row, "StrainId"]]], col] <- 
          samplingInfo[row, col]
      }
    }
  }
  
  return(tipInfo)
}

indexArray <- function(array){
  
  output <- list()
  for(i in 1:length(array)){
    output[[array[i]]] <- i
  }
  
  return(output)
}

removeBadgersFromTree <- function(tree){
  
  # Examine each of the tips on the phylogenetic tree
  for(tip in tree$tip.label){
    
    # Is the current tip for a badger isolate?
    if(grepl(pattern="WB", x=tip) == TRUE){
      tree <- drop.tip(tree, tip)
    }
  }
  
  return(tree)
}
