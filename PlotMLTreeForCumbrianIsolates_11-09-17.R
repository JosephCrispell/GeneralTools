# Load the ape package
library(ape)
library(geiger) # For the tips function
library(plotrix)

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Cumbria/"

###################################
# Get the Maximum Likelihood Tree #
###################################

# Read in the newick tree
file <- paste(path, "vcfFiles/",
              "mlTree_10-09-17.tree", sep="")
tree <- read.tree(file=file)

# Drop reference
tree <- drop.tip(tree, "Ref-1997")

# Convert Branch lengths to SNPs
fastaLength <- 2071
tree$edge.length <- tree$edge.length * fastaLength

##################################
# Get the Isolate Coverage Table #
##################################

file <- paste(path,"vcfFiles/",
              "IsolateVariantPositionCoverage_10-09-2017.txt", sep="")
table <- read.table(file, header=TRUE, stringsAsFactors=FALSE)

table$Isolate <- getIsolateIDFromFileNames(table$Isolate)

######################
# Get sampling dates #
######################

file <- paste(path, "17z_metadata_280717.csv", sep="")
samplingInfo <- read.table(file, header=TRUE, stringsAsFactors=FALSE,
                           sep=",")

samplingInfo$Isolate <- parseIsolateIDs(samplingInfo$Isolate)

# Make tip labels include sampling location and date
tips <- tree$tip.label
tree$tip.label <- changeTipLabels(tree$tip.label, samplingInfo)

##############################
# Plot the Phylogenetic Tree #
##############################

file <- paste(path, "vcfFiles/", "mlTree_10-09-17.pdf", sep="")
pdf(file, height=10, width=10)

# Set the margins
par(mfrow=c(1,1))
par(mar=c(0,0,0,0)) # Bottom, Left, Top, Right

plotType <- "fan" # "phylogram", "cladogram", "fan", "unrooted", "radial"

# Get each isolate's quality
isolateQuality <- getIsolateQuality(table)
factor <- 2

# Plot the phylogenetic tree
plot.phylo(tree, show.tip.label=TRUE, plotType,
           edge.color="dimgrey", edge.width=3,
           show.node.label=TRUE, label.offset=0.15,
           underscore=TRUE)

# Add node labels
nodelabels(node=1:length(tree$tip.label), 
           cex=defineTipSizeBySequencingQuality(tips, isolateQuality) * factor,
           pch=ifelse(tips %in% c("AF05293-17", "AF05220-17", "AF05712-17"),
                      21, 24),
           bg=ifelse(tips %in% c("AF05293-17", "AF05220-17", "AF05712-17"),
                     "red", "blue"), 
           col="dimgrey")

# Add Legends
text(x=5.25, y=-2.15, labels="Coverage:", col="black", cex=1)
addLegendForQuality("bottomright", 1, factor=factor)
legend("bottomleft", legend=c("Cow", "Badger"),
       pch=c(17, 16), cex=1, col=c("blue", "red"), 
       text.col=c("blue", "red"), bty='n')

# Add Scale bar
points(x=c(-0.5, 0.5), y=c(-4, -4), type="l", lwd=3)
text(x=0, y=-4.25, labels="~ 1 SNP", cex=1)


dev.off()

#############
# FUNCTIONS #
#############

changeTipLabels <- function(tipLabels, samplingInfo){
  output <- c()
  for(i in 1:length(tipLabels)){
    
    row <- which(samplingInfo$Isolate == tipLabels[i])
    
    output[i] <- paste(samplingInfo[row, "Location"],
                       samplingInfo[row, "Date"], sep="_")
  }
  
  return(output)
}

addLegendForQuality <- function(position, cex, factor){
  
  sizes <- seq(0, 1, 0.1)
  
  legend(position, legend=sizes, col="black", pch=24, bty='n',
         pt.cex=sizes * factor, cex=cex)
}

parseIsolateIDs <- function(names){
  
  isolates <- c()
  for(i in 1:length(names)){
    
    if(grepl(pattern="T-15-525", names[i]) == TRUE){
      isolates[i] <- names[i]
    }else{
      parts <- strsplit(names[i], split="-")[[1]]
      isolates[i] <- paste(parts[1], parts[3], "-", parts[4], sep="")
    }
  }
  
  return(isolates)
}

getIsolateIDFromFileNames <- function(fileNames){
  isolates <- c()
  for(i in 1:length(fileNames)){
    
    if(grepl(pattern="T-15-525", fileNames[i]) == TRUE){
      isolates[i] <- strsplit(fileNames[i], split="_")[[1]][1]
    }else{
      parts <- strsplit(strsplit(fileNames[i], split="_")[[1]][1], split="-")[[1]]
      isolates[i] <- paste(parts[1], parts[3], "-", parts[4], sep="")
    }
  }
  
  return(isolates)
}

getIsolateQuality <- function(table){
  isolateQuality <- list()
  for(i in 1:nrow(table)){
    isolateQuality[[table[i, "Isolate"]]] <- table[i, "VariantPositionCoverage"]
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
