#############
# Libraries #
#############

library(ape)
library(geiger)

###############################################
# Read in the inter isolate genetic distances #
###############################################

# Set the path
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_13-07-17/"

# Read in the table
file <- paste(path, "GeneticVsEpidemiologicalDistances/",
  "GeneticVsEpidemiologicalDistances_02-10-17.txt", sep="")
geneticVsEpi <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)

####################################
# Identify isolates in BASTA clade #
####################################

# Read in the newick tree
file <- paste(path, "vcfFiles/",  "mlTree_29-09-2017.tree", sep="")
isolatesInClade <- getIsolatesInClade(file, node=301)

######################################################
# Get distances only between isolates in BASTA clade #
######################################################

geneticVsEpi <- removeComparisonsBetweenIsolatesNotInClade(geneticVsEpi,
                                                           isolatesInClade)

##################################################################################
# Compare the distributions of BB, CC, and BC/CB inter-isolate genetic distances #
##################################################################################

file <- paste(path, "InterBASTAIsolateGeneticDistances_15-11-17.pdf", sep="")
pdf(file)

par(mar=c(5.1, 4.1, 4.1, 2.1))

boxplot(geneticVsEpi$GeneticDistance ~ geneticVsEpi$iSpeciesJSpecies,
        las=1, pch=20, outline=FALSE,
        ylab="Genetic Distance (N. SNPs)", col="white", border="white")

stripchart(geneticVsEpi$GeneticDistance ~ geneticVsEpi$iSpeciesJSpecies,
           vertical = TRUE, jitter=0.2,
           method = "jitter", add = TRUE, pch = 20,
           col=rgb(0.5,0.5,0.5, 0.1))

boxplot(geneticVsEpi$GeneticDistance ~ geneticVsEpi$iSpeciesJSpecies,
        las=1, pch=20, outline=FALSE, yaxt="n", xaxt="n",
        ylab="", border=rgb(0,0,0, 0.75), add=TRUE,
        col=rgb(0,0,0,0))

dev.off()


#############
# FUNCTIONS #
#############

removeComparisonsBetweenIsolatesNotInClade <- function(geneticVsEpi, isolatesInClade){
  
  keep <- c()
  
  for(row in 1:nrow(geneticVsEpi)){
    
    if(geneticVsEpi[row, "IsolateI"] %in% isolatesInClade &&
       geneticVsEpi[row, "IsolateJ"] %in% isolatesInClade){
      keep[length(keep) + 1] <- row
    }
  }
  
  return(geneticVsEpi[keep, ])
}

getIsolatesInClade <- function(treeFile, node){
  
  tree <- read.tree(file=treeFile)
  
  # Get a list of the isolates in the clade
  cladeTips <- tips(tree, node=node)
  
  return(cladeTips)
}
