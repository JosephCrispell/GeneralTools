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

########################################
# Define Clades by Threshold Distances #
########################################

# Note the nodes defining clades
nodesDefiningClades <- c(298, 276, 250, 487, 481, 460, 441, 291, 282, 261) # use nodelabels() to show node numbers

# Colour branches of clades
cladeColours <- c("red", "blue", "green", "cyan", "orange", "darkorchid4", "deeppink", "black", "brown", "darkolivegreen")
#cladeColours <- rep("red", length(nodesDefiningClades))
branchColours <- defineBranchColoursOfClades(tree, nodesDefiningClades, cladeColours, "lightgrey")

# Colour tips of clades
tipColours <- getTipColoursBasedOnClades(tree, nodesDefiningClades, cladeColours, "lightgrey")

# Define tip shapes by species
tipShapes <- defineTipShapesForSpecies(tree$tip.label, 24, 21)

##################################
# Get the Isolate XY coordinates #
##################################

path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/IsolateData/"

# Cattle Isolates
file <- paste(path, "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",")

# Badger Isolates
file <- paste(path, "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
badgerInfo <- read.table(file, header=TRUE, sep=",")

tipXYs <- matrix(nrow=length(tree$tip.label), ncol=2)
colnames(tipXYs) <- c("X", "Y")

for(tipIndex in 1:length(tree$tip.label)){
  
  # Cattle
  if(grepl(pattern="TB", x=tree$tip.label[tipIndex]) == TRUE){
    
    # Find index in table
    strainIndex <- which(cattleInfo$StrainId == tree$tip.label[tipIndex])
    tipXYs[tipIndex, 1] <- cattleInfo[strainIndex, "Mapx"]
    tipXYs[tipIndex, 2] <- cattleInfo[strainIndex, "Mapy"] 
  
  # Badgers
  }else if(grepl(pattern="WB", x=tree$tip.label[tipIndex]) == TRUE){

    # Find index in table
    strainIndex <- which(badgerInfo$WB_id == tree$tip.label[tipIndex])
    tipXYs[tipIndex, 1] <- badgerInfo[strainIndex, "SampledGrpX"]
    tipXYs[tipIndex, 2] <- badgerInfo[strainIndex, "SampledGrpY"] 
  }
}

####################################################
# Plot the Phylogenetic Tree and Isolate Locations #
####################################################

file <- paste(path, "mlTreeAndIsolateLocations_01-12-16.pdf", sep="")
pdf(file, height=7, width=14)

par(mfrow=c(1,2))

## Plot the phylogenetic tree
plot.phylo(tree, show.tip.label=FALSE, "fan",
           edge.color=branchColours, edge.width=3)

nodelabels(node=1:length(tree$tip.label),
           pch=tipShapes,
           col="dimgrey",
           cex=0.75,
           bg=tipColours)

points(x=c(100, 110), y=c(-90, -90), type="l", lwd=3)
text(x=105, y=-95, labels="10 SNPs", cex=0.75)
mtext("A", side=3, at=-116)

legend(x=-120, y=-85, legend=c("Cow", "Badger"),
       pch=c(17, 16), cex=0.65, bty='n')

## Plot the Isolate Locations

# Note Mansion Location
mansionX <- 380909
mansionY <- 201377

# Set an expansion
expand <- 12500

# Don't plot unassigned
tipColours[tipColours == "lightgrey"] <- NA

# Create empty plot
plot(1, type="n", yaxt="n", xaxt="n",
     xlim=c(mansionX - expand, mansionX + expand),
     ylim=c(mansionY - expand, mansionY + expand),
     xlab=paste((expand * 2) / 1000, "km"), ylab=paste((expand * 2) / 1000, "km"),
     main="Isolate Locations")

legend("topleft", legend=c("Cow", "Badger", "Woodchester Mansion"),
       pch=c(17, 16, 15), cex=0.65, bty='n')

# Add point for Woodchester Mansion
points(x=mansionX, y=mansionY, pch=15, col="black")

# Add polygons around clade points
addPolygonsForClades(cladeColours, tipXYs, tipColours, 0.35)

# Set transparency
alpha = 0.8

# Add the isolates
for(tipIndex in 1:length(tree$tip.label)){
  
  if(is.na(tipColours[tipIndex]) == TRUE){
    next
  }
  
  points(x=tipXYs[tipIndex, "X"], y=tipXYs[tipIndex, "Y"], 
         pch=tipShapes[tipIndex], cex=0.75,
         col=setColourTransparency(tipColours[tipIndex], alpha), 
         bg=setColourTransparency(tipColours[tipIndex], alpha))
}
mtext("B", side=3, at=368000)
dev.off()

#############
# FUNCTIONS #
#############

addPolygonsForClades <- function(cladeColours, tipXYs, tipColours, alpha){
  
  cladePoints <- as.data.frame(tipXYs)
  cladePoints$Colours <- tipColours
  
  # Remove unassigned
  cladePoints <- cladePoints[is.na(cladePoints$Colours) == FALSE, ]
  
  # Remove no locations
  cladePoints <- cladePoints[is.na(cladePoints$X) == FALSE, ]
  
  for(colour in cladeColours){
    subset <- cladePoints[cladePoints$Colours == colour, ]
    addPolygon(xValues = subset$X, yValues = subset$Y, 
               borderColour = setColourTransparency(colour, alpha))
  }
}

addPolygon <- function(xValues, yValues, borderColour){
  hullPoints <- chull(xValues, yValues)
  hullPoints <- c(hullPoints, hullPoints[1])
  polygon(xValues[hullPoints], yValues[hullPoints], col = NA, border = borderColour)
}

setColourTransparency <- function(color, alpha) {

  # Set the percentage transparency
  percent <- (1-alpha)*100
  
  # Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  # Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100)

  return(t.col)
}

getTipColoursBasedOnClades <- function(tree, nodesDefiningClades, cladeColours, defaultColour){
  tipColours <- rep(defaultColour, length(tree$tip.label))
  tipsInClades <- getTipsInClades(tree, nodesDefiningClades)
  for(tipIndex in 1:length(tree$tip.label)){
    
    for(nodeIndex in 1:length(nodesDefiningClades)){
      
      if(tree$tip.label[tipIndex] %in% tipsInClades[[nodesDefiningClades[nodeIndex]]] == TRUE){
        tipColours[tipIndex] <- cladeColours[nodeIndex]
        break
      }
    }
  }
  
  return(tipColours)
}

getTipsInClades <- function(tree, nodesDefiningClades){
  tipsInClades <- list()
  for(node in nodesDefiningClades){
    tipsInClades[[node]] <- tips(tree, node)
  }
  return(tipsInClades)
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


defineBranchColoursOfClades <- function(tree, nodesDefiningClades,
                                        CladeColours, defaultColour){
  branchColours <- rep(defaultColour, dim(tree$edge)[1])
  for(i in 1:length(cladeColours)){
    clade <- tips(tree, node=nodesDefiningClades[i])
    branchesInClades <- which.edge(tree, clade)
    branchColours[branchesInClades] <- cladeColours[i]
  }
  return(branchColours)
}