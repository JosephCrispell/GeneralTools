# Load the ape package
library(ape)
library(geiger) # For the tips function

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/"

################################
# Read in sampling information #
################################

# Read in the badger isolate metadata
fileName <- paste(path, "IsolateData/", "BadgerInfo_08-04-15_LatLongs_XY.csv", sep="")
badgerInfo <- read.table(fileName, header=TRUE, stringsAsFactors=FALSE, sep=",")

# Read in the cattle isolate metadata
file <- paste(path, "IsolateData/", 
              "CattleIsolateInfo_LatLongs_plusID_outbreakSize_Coverage_AddedTB1453-TB1456.csv", sep="")
cattleInfo <- read.table(file, header=TRUE, sep=",", stringsAsFactors=FALSE)

##################################
# Get the Isolate Coverage Table #
##################################

file <- paste(path,"allVCFs-IncludingPoor/vcfFiles/",
              "isolateGenomeCoverageSummary_28-10-16.txt", sep="")
coverageInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

coverageInfo$IsolateID <- getIsolateID(coverageInfo$IsolateID)

#############################################
# Get the isolate variant position coverage #
#############################################

file <- paste(path,"allVCFs-IncludingPoor/vcfFiles/",
              "IsolateVariantPositionCoverage_.txt", sep="")
vpCoverage <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

vpCoverage$Isolate <- getIsolateID(vpCoverage$Isolate)

#######################################
# Open the tip to root distances file #
#######################################

# Read in the file
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "DistancesToRoot_mlTree_28-10-16.txt", sep="")
tipToRootDistances <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

# Add sampling times
tipToRootDistances$SamplingTime <- getSamplingTimes(tipToRootDistances$IsolateID,
                                                    badgerInfo, cattleInfo)

# Add genome coverage
tipToRootDistances$GenomeCoverage <- getIsolateGenomeCoverage(
  tipToRootDistances$IsolateID, coverageInfo)

# Add vp coverage
tipToRootDistances$VariantPositionCoverage <- getIsolateVPCoverage(
  tipToRootDistances$IsolateID, vpCoverage)

###################################
# Get the maximum likelihood tree #
###################################

# Read in the newick tree
file <- paste(path, "allVCFs-IncludingPoor/vcfFiles/",
              "mlTree_Prox-10_plusRef_rmResequenced_SNPCov-0.1_28-10-16.tree", sep="")
tree <- read.tree(file=file)

# Convert Branch lengths to SNPs
fastaLength <- 9464
tree$edge.length <- tree$edge.length * fastaLength

###############################
# Select clade used for BASTA #
###############################

# Set margins
par(mar=c(0,0,0,0)) # Bottom, Left, Top, Right (Defaults: 5.1,4.1,4.1,2.1)

# Get the BASTA clade
node <- 289
bastaClade <- extract.clade(tree, node)

# Plot the BASTA clade within full tree
branchColours <- rep("black", dim(tree$edge)[1])
branchColours[which.edge(tree, tips(tree, node=node))] <- "red"
plot.phylo(tree, show.tip.label=FALSE, "fan",
           edge.color=branchColours, edge.width=3)

# Plot BASTA clade alone
plot.phylo(bastaClade, show.tip.label=FALSE, "fan", edge.width=3)

###########################
# Start looking at clades #
###########################

# Create a list of potential heterozygous isolates
heterozygous <- c("TB1775", "TB1810", "TB1811", "TB1391", "TB1415")
tipColours <- rep("black", length(bastaClade$tip.label))
tipColours[bastaClade$tip.label %in% heterozygous] <- "red"
plot.phylo(bastaClade, show.tip.label=TRUE, "fan", edge.width=3,
           tip.color=tipColours, label.offset=0.1, cex=0.75)

tipColours <- rep("black", length(tree$tip.label))
tipColours[tree$tip.label %in% heterozygous] <- "red"
plot.phylo(tree, show.tip.label=TRUE, "fan", edge.width=3,
           tip.color=tipColours, label.offset=3, cex=0.6)


# Create a list of flagged isolates on tree
flagged <- c("TB1501", "TB1503", "TB1505", "TB1776", "TB1783", "TB1390", "TB1393",
             "TB1402", "TB1404", "TB1405", "TB1408", "TB1448", "TB1452", "TB1458",
             "TB1463", "TB1485")
tipColours <- rep("black", length(bastaClade$tip.label))
tipColours[bastaClade$tip.label %in% flagged] <- "red"
plot.phylo(bastaClade, show.tip.label=TRUE, "fan", edge.width=3,
           tip.color=tipColours, label.offset=0.1, cex=0.75)

tipColours <- rep("black", length(tree$tip.label))
tipColours[tree$tip.label %in% flagged] <- "red"
plot.phylo(tree, show.tip.label=TRUE, "fan", edge.width=3,
           tip.color=tipColours, label.offset=3, cex=0.6)

# Colour tips by time
samplingTimes <- getSamplingTimes(bastaClade$tip.label, badgerInfo, cattleInfo)

rbPal <- colorRampPalette(c('blue','red'))
tipColours <- rbPal(10)[as.numeric(cut(samplingTimes,breaks = 10))]

plot.phylo(bastaClade, show.tip.label=TRUE, "fan", edge.width=3,
           tip.color=tipColours, label.offset=0.15, cex=0.75,
           show.node.label=TRUE)

nodelabels(node=1:length(bastaClade$tip.label), 
           cex=1, pch=16,
           col=tipColours)

times <- seq(min(samplingTimes), max(samplingTimes), 1000)
timeColours <- rbPal(10)[as.numeric(cut(times,breaks = 10))]

legend(x=15, y=15, legend=as.Date(times, origin="1970-01-01"),
       col=timeColours, pch=16, bty='n')

# Define clades within tree
nodesDefiningClades <- c(229, 204, 365, 345, 317, 196) # use nodelabels(frame="none", col="red", font=2) to show node numbers

# Colour branches of clades
cladeColours <- c("red", "blue", "green", "cyan", "orange", "darkorchid4")
branchColours <- defineBranchColoursOfClades(bastaClade, nodesDefiningClades,
                                             cladeColours, "lightgrey")
plot.phylo(bastaClade, show.tip.label=TRUE, "fan", edge.width=3,
           label.offset=0.15, cex=0.75, edge.color=branchColours)

# Note the clades of each tip
tipClades <- noteCladesOfTips(bastaClade, nodesDefiningClades)

# file <- "C:/Users/Joseph Crisp/Desktop/tree.pdf"
# pdf(file, height=20, width=20)
# plot.phylo(bastaClade, show.tip.label=TRUE, "fan", edge.width=3,
#            label.offset=0.15, cex=0.75, show.node.label=TRUE)
# nodelabels(frame="none", col="red", font=2)
# dev.off()

# Add clade colours to root-to-tip table
tipToRootDistances$CladeColours <- getIsolateCladeColours(tipToRootDistances$IsolateID,
                                                          tipClades, cladeColours,
                                                          "lightgrey")

# Select only the BASTA clade isolates
tipToRootDistances <- tipToRootDistances[
  tipToRootDistances$IsolateID %in% bastaClade$tip.label, ]

# Plot the model output
par(mar=c(5.1,4.1,4.1,2.1))
plot(tipToRootDistances$SamplingTime, tipToRootDistances$PatristicDistanceToRoot,
     las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=tipToRootDistances$CladeColours, cex=2,
     pch=ifelse(grepl(x=tipToRootDistances$IsolateID, pattern="WB") == TRUE,
                19, 17))

for(colour in cladeColours){
  abline(lm(PatristicDistanceToRoot~SamplingTime, 
            data=tipToRootDistances[tipToRootDistances$CladeColours == colour, ]),
         col=colour)
}

# Plot the root-to-tip correlation and colour by genome coverage
rbPal <- colorRampPalette(c('blue','red'))
colours <- rbPal(10)[as.numeric(cut(tipToRootDistances$GenomeCoverage,breaks = 10))]
plot(tipToRootDistances$SamplingTime, tipToRootDistances$PatristicDistanceToRoot,
     las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=colours, cex=2,
     pch=ifelse(grepl(x=tipToRootDistances$IsolateID, pattern="WB") == TRUE,
                19, 17))

# Plot the root-to-tip correlation and colour by variant position coverage
plot(tipToRootDistances$GenomeCoverage, tipToRootDistances$VariantPositionCoverage,
     las=1, ylab="Variant Position Coverage", xlab="Genome Coverage",
     xlim=c(0,1), ylim=c(0,1), col=rgb(0,0,0, 0.5),
     pch=ifelse(grepl(x=tipToRootDistances$IsolateID, pattern="WB") == TRUE,
                                                                     19, 17))

rbPal <- colorRampPalette(c('blue','red'))
colours <- rbPal(10)[as.numeric(cut(tipToRootDistances$VariantPositionCoverage,breaks = 10))]
plot(tipToRootDistances$SamplingTime, tipToRootDistances$PatristicDistanceToRoot,
     las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=colours, cex=2,
     pch=ifelse(grepl(x=tipToRootDistances$IsolateID, pattern="WB") == TRUE,
                19, 17))

# Account for coverage in root-to-tip correlation
rbPal <- colorRampPalette(c('blue','red'))
colours <- rbPal(10)[as.numeric(cut(tipToRootDistances$VariantPositionCoverage,breaks = 10))]
plot(tipToRootDistances$SamplingTime, 
     tipToRootDistances$PatristicDistanceToRoot * (1/tipToRootDistances$VariantPositionCoverage),
     las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=colours, cex=2,
     pch=ifelse(grepl(x=tipToRootDistances$IsolateID, pattern="WB") == TRUE,
                19, 17))

linearModel <- lm(tipToRootDistances$PatristicDistanceToRoot * (1/tipToRootDistances$VariantPositionCoverage) ~ tipToRootDistances$SamplingTime)
summary <- summary(linearModel)
abline(linearModel, col="black", lwd=2)

legend("topright", 
       legend=c(paste("R^2 = ", round(summary$adj.r.squared, 2)),
                paste("p-value = ", round(anova(linearModel)$Pr[[1]], 2))),
       bty='n')

# Look at removing poor coverage isolates
subset <- tipToRootDistances[tipToRootDistances$VariantPositionCoverage > 0.8, ]

rbPal <- colorRampPalette(c('blue','red'))
colours <- rbPal(10)[as.numeric(cut(tipToRootDistances$VariantPositionCoverage,breaks = 10))]
plot(subset$SamplingTime, subset$PatristicDistanceToRoot,
     las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=colours, cex=2,
     pch=ifelse(grepl(x=subset$IsolateID, pattern="WB") == TRUE,
                19, 17))

linearModel <- lm(subset$PatristicDistanceToRoot ~ subset$SamplingTime)
summary <- summary(linearModel)
abline(linearModel, col="black", lwd=2)

legend("topright", 
       legend=c(paste("R^2 = ", round(summary$adj.r.squared, 2)),
                paste("p-value = ", round(anova(linearModel)$Pr[[1]], 2))),
       bty='n')

rbPal <- colorRampPalette(c('blue','red'))
colours <- rbPal(10)[as.numeric(cut(tipToRootDistances$VariantPositionCoverage,breaks = 10))]
plot(subset$SamplingTime, 
     subset$PatristicDistanceToRoot * (1/subset$VariantPositionCoverage),
     las=1, ylab="Root-to-tip Distance", xlab="Sampling Date",
     main="", mgp=c(3,0.5,0), col=colours, cex=2,
     pch=ifelse(grepl(x=subset$IsolateID, pattern="WB") == TRUE,
                19, 17))

linearModel <- lm(subset$PatristicDistanceToRoot * (1/subset$VariantPositionCoverage) ~ subset$SamplingTime)
summary <- summary(linearModel)
abline(linearModel, col="black", lwd=2)

legend("topright", 
       legend=c(paste("R^2 =", round(summary$adj.r.squared, 2)),
                paste("p-value =", round(anova(linearModel)$Pr[[1]], 2)),
                "Coverage threshold = 0.8"),
       bty='n')

linearModel <- lm(subset$PatristicDistanceToRoot * (1/subset$VariantPositionCoverage) ~ 
                    subset$SamplingTime + subset$VariantPositionCoverage)
summary <- summary(linearModel)


#############
# FUNCTIONS #
#############

getIsolateVPCoverage <- function(isolateIDs, vpCoverage){
  
  isolateCoverage <- c()
  for(i in 1:length(isolateIDs)){
    
    row <- which(vpCoverage$Isolate == isolateIDs[i])
    if(length(row) != 0){
      isolateCoverage[i] <- vpCoverage[row, "VariantPositionCoverage"]
    }
  }
  
  return(isolateCoverage)
}

getIsolateGenomeCoverage <- function(isolateIDs, coverageInfo){
  
  isolateCoverage <- c()
  for(i in 1:length(isolateIDs)){
    
    row <- which(coverageInfo$IsolateID == isolateIDs[i])
    if(length(row) != 0){
      isolateCoverage[i] <- coverageInfo[row, "PercentageCoverage"]
    }
  }
  
  return(isolateCoverage)
}

getIsolateID <- function(fileNames){
  ids <- c()
  for(i in 1:length(fileNames)){
    ids[i] <- strsplit(fileNames[i], split="_")[[1]][1]
  }
  
  return(ids)
}

getIsolateCladeColours <- function(isolateIDs, tipClades, cladeColours, default){
  
  isolateCladeColours <- rep(default, length(isolateIDs))
  
  for(i in 1:length(isolateIDs)){
    
    if(is.null(tipClades[[isolateIDs[i]]]) == FALSE){
      isolateCladeColours[i] <- cladeColours[tipClades[[isolateIDs[i]]]]
    }
  }
  
  return(isolateCladeColours)
}

noteCladesOfTips <- function(tree, nodesDefiningClades){

  tipClades <- list()
    
  for(i in 1:length(nodesDefiningClades)){
    
    for(tip in tips(tree, nodesDefiningClades[i])){
      tipClades[[tip]] <- i
    }
  }
  
  return(tipClades)
}

getSamplingTimes <- function(tips, badgerInfo, cattleInfo){
  
  # Initialise an array to store the sampling times associated with each tip
  samplingTimes <- c()
  
  # Note the format of the dates
  dateFormat <- "%d/%m/%Y"
  
  # Examine each tip
  for(i in 1:length(tips)){
    
    # Check if badger or cow isolate
    if(grepl(x=tips[i], pattern="WB") == TRUE){
      
      row <- which(badgerInfo$WB_id == tips[i])
      
      samplingTimes[i] <- as.Date(badgerInfo[row, "date"], dateFormat)
      
    }else{
      
      row <- which(cattleInfo$StrainId == tips[i])
      
      if(length(row) != 0){
        samplingTimes[i] <- as.Date(cattleInfo[row, "DateCultured"], dateFormat)
      }else{
        samplingTimes[i] <- NA
      }
    }
  }
  
  return(samplingTimes)
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
