# Load the ape package
library(ape)
library(geiger) # For the tips function

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/allVCFs-IncludingPoor/vcfFiles/"

###################################
# Get the Maximum Likelihood Tree #
###################################

# Read in the newick tree
file <- paste(path, "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Convert Branch lengths to SNPs
fastaLength <- 9464
tree$edge.length <- tree$edge.length * fastaLength

##################################
# Get the Isolate Coverage Table #
##################################

file <- paste(path, "isolateGenomeCoverageSummary_28-10-16.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

table$IsolateID <- getIsolateIDFromFileNames(table$IsolateID)

##############################
# Plot the Phylogenetic Tree #
##############################

file <- paste(path, "mlTree_Prox-10_plusRef_rmReseq_SNPCov-0.1_Labelled-Clades-Species-Quality_28-10-16.pdf")
pdf(file, height=10, width=10)

plotType <- "fan" # "phylogram", "cladogram", "fan", "unrooted", "radial"

# Plot initial tree to find nodes defining clades
#plot.phylo(tree, plotType)
#nodelabels()

# Define branch colours by clade
#nodesDefiningClades <- c(515, 318, 521, 405, 433) # use nodelabels() to show node numbers
nodesDefiningClades <- c(291, 440, 460, 412, 324) # use nodelabels() to show node numbers
cladeColours <- c("cyan", "pink", "green", "orange", "purple")
branchColours <- defineBranchColoursOfClades(tree, nodesDefiningClades, cladeColours)

# Get each isolate's quality
isolateQuality <- getIsolateQuality(table)

# Plot the phylogenetic tree
plot.phylo(tree, show.tip.label=FALSE, plotType,
           edge.color=branchColours, edge.width=3,
           show.node.label=TRUE)

# Add node labels
nodelabels(node=1:length(tree$tip.label), 
           cex=defineTipSizeBySequencingQuality(tree$tip.label, isolateQuality),
           pch=defineTipShapesForSpecies(tree$tip.label, 17, 16),
           col=defineTipColourBySpecies(tree$tip.label, "blue", "red"))

# Add Legends
text(x=132, y=-68, labels="Coverage:", col="black", cex=0.7)
addLegendForQuality("bottomright", 0.8)
text(x=-109, y=-108, labels="Species:", col="black", cex=0.7)
legend("bottomleft", legend=c("COW", "BADGER"),
       pch=c(17, 16), cex=0.65, col=c("blue", "red"), 
       text.col=c("blue", "red"), bty='n')

# Add Scale bar
points(x=c(-20, 30), y=c(-114, -114), type="l", lwd=3)
text(x=5, y=-118, labels="50 SNPs", cex=0.8)

# Add Clade labels
text(x=52, y=102, labels="0", col=cladeColours[1], cex=2)
text(x=55, y=-100, labels="1", col=cladeColours[2], cex=2)
text(x=95, y=-75, labels="2", col=cladeColours[3], cex=2)
text(x=-45, y=-105, labels="3", col=cladeColours[4], cex=2)
text(x=-102, y=50, labels="4", col=cladeColours[5], cex=2)


dev.off()

#############
# FUNCTIONS #
#############

defineTipColourBySpecies <- function(tipLabels, cow, badger){
  colours <- c()
  for(i in 1:length(tipLabels)){
    
    if(grepl(pattern="TB", x=tipLabels[i]) == TRUE){
      colours[i] <- cow
    }else if(grepl(pattern="WB", x=tipLabels[i]) == TRUE){
      colours[i] <- badger
    }else{
      colours[i] <- "black"
    }
  }
  
  return(colours)
}

addLegendForQuality <- function(position, cex){

  sizes <- seq(0.1, 1, 0.1)
  
  legend(position, legend=sizes, col="black", pch=16, bty='n',
         pt.cex=sizes, cex=cex)
}

getIsolateQuality <- function(table){
  isolateQuality <- list()
  for(i in 1:nrow(table)){
    isolateQuality[[table[i, "IsolateID"]]] <- table[i, "PercentageCoverage"]
  }
  
  return(isolateQuality)
}

defineTipSizeBySequencingQuality <- function(tipLabels, isolateQuality){
  
  tipQuality <- c()
  for(i in 1:length(tipLabels)){
    if(tipLabels[i] != "Ref-1997"){
      
      tipQuality[i] <- isolateQuality[[tipLabels[i]]]
    }else{
      tipQuality[i] <- 1
    }
  }
  
  return(tipQuality)
}

defineBranchColoursOfClades <- function(tree, nodesDefiningClades,
                                        CladeColours){
  branchColours <- rep("black", dim(tree$edge)[1])
  for(i in 1:length(cladeColours)){
    clade <- tips(tree, node=nodesDefiningClades[i])
    branchesInClades <- which.edge(tree, clade)
    branchColours[branchesInClades] <- cladeColours[i]
  }
  return(branchColours)
}

getUpperBoundsFromCuts <- function(cuts){
  
  bounds <- c()
  for(i in 1:length(cuts)){
    bounds[i] <- as.numeric(strsplit(strsplit(cuts[i], ",")[[1]][2], "]")[[1]][1])
  }
  
  return(bounds)
}

addLegend <- function(position, colours, nBreaks, cex){
  
  colourPalette <- colorRampPalette(colours)
  
  cuts <- levels(cut(table$PercentageCoverage, breaks=nBreaks))
  
  bounds <- getUpperBoundsFromCuts(cuts)
  
  bounds <- round(bounds, 2)
  
  legend(position, legend=bounds, col=colourPalette(nBreaks), pch=20, bty='n',
         cex=cex)
  
}

getIsolateIDFromFileNames <- function(fileNames){
  isolates <- c()
  for(i in 1:length(fileNames)){
    isolates[i] <- strsplit(fileNames[i], split="_")[[1]][1]
  }
  
  return(isolates)
}

assignIsolatesContinuousColoursByCoverage <- function(table, colours, nBreaks){
  
  colourPalette <- colorRampPalette(colours)
  coloursPerRow <- colourPalette(nBreaks)[
    as.numeric(cut(table$PercentageCoverage, breaks=nBreaks))]
  
  isolateColours <- list()
  for(i in 1:nrow(table)){
    isolateColours[[table[i, 1]]] <- coloursPerRow[i]
  }
  
  return(isolateColours)
}

returnTipColoursForIsolates <- function(tipLabels, assignedColours){
  
  tipColours <- c()
  for(i in 1:length(tipLabels)){
    if(tipLabels[i] != "Ref-1997"){
      
      tipColours[i] <- assignedColours[[tipLabels[i]]]
    }else{
      tipColours[i] <- "black"
    }
  }
  
  return(tipColours)
}

defineTipShapesForSpecies <- function(tipLabels, cow, badger){
  
  shapes <- c()
  for(i in 1:length(tipLabels)){
    
    if(grepl(pattern="TB", x=tipLabels[i]) == TRUE){
      shapes[i] <- cow
    }else if(grepl(pattern="WB", x=tipLabels[i]) == TRUE){
      shapes[i] <- badger
    }else{
      shapes[i] <- cow
    }
  }
  
  return(shapes)
}
