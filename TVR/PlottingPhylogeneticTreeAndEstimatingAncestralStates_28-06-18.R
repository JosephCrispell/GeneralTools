#### PACKAGES ####
library(ape)
library(phangorn)
library(phytools)

#### FUNCTIONS ####

countTransitions <- function(tree, tipStates, ancestralStateProbs, probThreshold=0.5){
  
  # CURRENTLY TOO SIMPLISTIC! FORGOT THAT NEED TO CONSIDER CHERRIES!!!
  # A more conservative transition count method - recommended by Nicola De Maio
  # 	- Assumes the ancestor represents one of the tips in the past
  # 
  # For a two state problem: A & B
  #  ----A				----A				----A				  ----B	
  # |		  1 AA	 |		0 AA	 |		  1 AA	 |		 0 AA
  # A		  0 AB	 A		1 AB	 A----B	1 AB	 A		 2 AB
  # |				     |				   |				     |
  # ----A			   ----B			 ----A				 ----B
  
  # Nodes in phylogenetic tree are numbered:
  # 1:Ntip(tree) (1:Ntip(tree) + tree$Nnode)
  # 
  # Edges on the phylogenetic tree are under tree$edge
  # A matrix of two columns: FROM, TO - based upon the node indices defined above
  #
  # Ancestral state probs - probability of each state for each INTERNAL node
  #
  # Probability threshold is the threshold to define the state assigned to each internal node
  
  # Note the state probabilities for the tips
  tipStateProbs <- matrix(0, nrow=length(tipStates), ncol=2)
  tipStateProbs[tipStates == "Bovine", 2] <- 1
  tipStateProbs[tipStates == "Badger", 1] <- 1
  
  # Add the state probabilities for the tips onto the ancestral state probs for ease
  stateProbs <- rbind(tipStateProbs, ancestralStateProbs)
  
  # Initialise a list to store the states associated with each internal nodes daughter nodes
  nodes <- list()
  
  # Examine each edge in the phylogeny
  for(row in 1:nrow(tree$edge)){
    
    # Note the node index of from node
    fromNode <- as.character(tree$edge[row, 1])
    
    # Get the state for the from node
    fromState <- which(stateProbs[as.numeric(fromNode), ] > probThreshold)
    if(length(fromState) == 0){
      fromState <- -1
    }
    
    # Get the state for the to node
    toState <- which(stateProbs[tree$edge[row, 2], ] > probThreshold)
    if(length(toState) == 0){
      toState <- -1
    }
    
    # Check if we have encountered this from node before
    if(is.null(nodes[[fromNode]]) == FALSE){
      
      nodes[[fromNode]] <- c(nodes[[fromNode]], toState)
      
    # Create a record for the current node if we haven't encountered it before
    }else{
      nodes[[fromNode]] <- c(fromState, toState)
    }
  }
  
  # Initialise a matrix to count the transitions
  transitionCounts <- matrix(0, nrow=2, ncol=2)
  colnames(transitionCounts) <- colnames(stateProbs)
  rownames(transitionCounts) <- colnames(stateProbs)
  
  # Examine each of the node
  for(node in names(nodes)){
    
    # Initialise a variable to count the number of branches beginning and ending with the same state
    nSameStateBranches <- 0
    
    # Skip node if no state available
    if(nodes[[node]][1] == -1){
      next
    }
    
    # Examine each of the daughter states
    for(i in 2:length(nodes[[node]])){
      
      # Check if state of daughter is the same as that of parent
      if(nodes[[node]][i] != -1 && nodes[[node]][i] == nodes[[node]][1]){
        nSameStateBranches <- nSameStateBranches + 1
        
      # Found branch beginning and ending with different states - count the inter species transmission event
      }else if(nodes[[node]][i] != -1){
        transitionCounts[nodes[[node]][1], nodes[[node]][i]] <- transitionCounts[nodes[[node]][1], nodes[[node]][i]] + 1
      }
    }
    
    # Count the number of within species transmission events - conservatively assumes parent node is same animal as one of daughters
    if(nSameStateBranches > 1){
      for(i in seq_len(nSameStateBranches - 1)){
        transitionCounts[nodes[[node]][1], nodes[[node]][1]] <- transitionCounts[nodes[[node]][1], nodes[[node]][1]] + 1
      }
    }
  }
  
  return(transitionCounts)
}

#### READ FASTA ####

# Read in the FASTA file
sequences<- read.dna("/home/josephcrispell/Desktop/Research/TVR/263_single_isolate_Cattle&Badgers_all_years.fasta",
                            format = "fasta")

# Note the number of sites in the FASTA file
nSitesInFasta <- length(as.list(sequences)[[1]])

#### CHOOSE SUBSTITUTION MODEL ####

# Convert the sequences from a DNAbin format into a PhyDat format
sequencesPhyDat <- phyDat(sequences, type = "DNA", levels = NULL)

# Run model testing to select appropriate substitution model
modelTestResults <- modelTest(sequencesPhyDat, model = c("JC", "HKY", "GTR"))

# Get best model
cat(paste("Best substitution model:", 
          modelTestResults$Model[which.min(modelTestResults$AIC)], "\n"))


#### BUILD PHYLOGENY ####

# Build the distance matrix
distanceMatrix <- dist.dna(sequences, model="JC69")

# Build a neighbour joining tree - an initial tree
njTree <- nj(distanceMatrix)

# Compute likelihood of the initial Neighbour Joining tree given sequences
likelihoodObject <- pml(njTree, sequencesPhyDat)

# Set the controls for the Maximum Likelihood algorithm
controls <- pml.control(maxit=100000, trace=0)

# Run maximum likelihood
fittingOutput <- optim.pml(likelihoodObject,
                           optNni = TRUE, # Optimise topology
                           optInv = TRUE, # Optimise proportion of variable sites
                           model = "GTR", # Substitution model
                           rearrangement="NNI", # Nearest Neighbour Interchanges
                           control=controls)

# Bootstrap the result of maximum likelihood
bootstrapResults <- bootstrap.pml(
  fittingOutput, # Use Maximium Likelihood settings on bootstrapped sequences
  bs = 100, # Number times to bootstrap sequences
  optNni = TRUE, # Use Nearest Neighbour Interchanges in tree building
  jumble=TRUE) # Jumble bootstrapped sequences before building trees

# Get phylogenetic tree with bootstrap values - plotBS will return a tree that can be used
tree <- plotBS(
  fittingOutput$tree,
  bootstrapResults,
  p = 50, # Plot bootstrap values if node in >=50 bootstrap trees
  type="phylogram", cex=0.1) # Type of phylogenetic tree shape to plot

### To ensure you can get a SNP difference scale bar on subsequent phylogenies.
### You need to alter the branch lengths to SNPs and not proportions.
tree$edge.length<-tree$edge.length*nSitesInFasta

###If you want to write the tree to the Newick format that can be read by Figtree, do the following.
write.tree(phy=tree, file="TVR phylogeny Feb19 - n=303 single isolates all years")


#### ANCESTRAL STATE RECONSTRUCTION ####

# Remove any non-dichotomous branching nodes (nodes with more than two daughters)
rootedTree <- multi2di(tree, random=TRUE)

# Root the tree
rootedTree <- root(rootedTree, 
                   outgroup="Reference", # Set the isolate to root the tree on - can also use node
                   resolve.root=TRUE) # Ensures new root will be bifurcating

# Make branches with lengths of 0 very very small
rootedTree$edge.length[rootedTree$edge.length <= 0] <- 0.000001

# Get tip states from their names
tipStates <- c()
for(index in 1:length(rootedTree$tip.label)){
  tipStates[index] <- strsplit(rootedTree$tip.label[index], split="_")[[1]][2]
}


### The lead state is the "Reference" genome - it has been entered as an NA because the scripts above were only pulling out the label badger or bovine.
### So, relabel this using the following script.
tipStates[1]<-"Bovine"

# Estimate ancestral states under a ER model
ancestralStateFitting <- ace(tipStates, rootedTree, 
                             model="ARD", # "ER" Equal rates; "ARD" All rates different
                             type="discrete") # Discrete states


#### PLOT TREE WITH ANCESTRAL STATES ####

# Get and set the plotting margins: c(bottom, left, top, right)
currentMar <- par()$mar
par(mar=c(0, 0, 0, 0))

# Assign colours to the states
stateColours <- setNames(c("red", "blue"), c("Badger", "Bovine"))

### Plot phylogenetic tree without ancestral states
plot.phylo(rootedTree, cex=0.1)

# Add tip shapes and colours - for the states
tiplabels(pie=to.matrix(tipStates,c("Badger", "Bovine")), # Provide the state of each tip
          piecol=stateColours, # Provide the colour assigned to each state
          cex=0.1) # Set the size of the tip circle

### Plot phylogenetic tree with ancestral states
plot.phylo(rootedTree, type="phylogram", # Set the shape of the tree to be a fan
           show.tip.label=FALSE) # Don't plot the tip labels

### Add scale bar
add.scale.bar(30, 0, length=2, lcol="red", lwd=1, cex=0.5)


# Add pie charts at each internal node
nodelabels(node= 1:rootedTree$Nnode+Ntip(rootedTree), # Provide an array with all the internal node IDs (nodes numbered from nTips onwards...)
           pie=ancestralStateFitting$lik.anc, # Assign pie charts to each node - based on the outputs from ancestral state estimation
           piecol=stateColours, # Set the colours to be used in the pie charts
           cex=0.1) # Set the size of the pie charts

# Add circles to indicate the states of the tips
tiplabels(pie=to.matrix(tipStates,c("Badger", "Bovine")), # Provide the state of each tip
          piecol=stateColours, # Provide the colour assigned to each state
          cex=0.1) # Set the size of the tip circle

# Add a legend
legend("bottomright",
       legend=c("Badger", "Bovine"), # Name the labels used in the legend
       pch=19, # Set the shape of the points used in legend
       col=c("red", "blue"), # Set the colour of each label point
       bty="n") # Remove box surrounding legend

#### COUNT TRANSITIONS ####

# Count the transitions based upon rooted tree with estimated ancestral states
transitionCounts <- countTransitions(tree=rootedTree, tipStates=tipStates,
                                     ancestralStateProbs=ancestralStateFitting$lik.anc,
                                     probThreshold=0.5)

#### Plot the transition counts ####

# Reset the plotting margins
par(mar=currentMar)

# Create a bar plot to show the estimate transition counts
barplot(c(transitionCounts[1,1], transitionCounts[1,2], transitionCounts[2,1], transitionCounts[2,2]),
     ylab="Number of estimate transitions", names=c("BB", "BC", "CB", "CC"), las=1)
