#################
# Load Packages #
#################

library(randomForest)

##################################################
# Read table - response versus predictor metrics #
##################################################

# Read in response vs predictor table
file <- "geneticVsPredictors_17-10-2017.txt"
table <- read.table(file, 
                    header=TRUE, # Input file has header
                    sep="\t") # coloumn separator in file is TAB

###########################
# Fit Random Forest model #
###########################

# Define the number of trees to be used in Random Forest
# You are basically looking for a number of trees that
# isn't so many that it takes ages but isn't too few 
# that the results from multiple runs aren't consistent
nTrees <- 1000

# Run the Random Forest tuning algorithm to find the optimal mtry value
# Mtry is the parameter that defines how many predictor variables to
# select from to build each node of each decision tree
# Random Forest has a built in tuning program to help you choose the
# right value for your data
tuneOutput <- tuneRF(table[-which(colnames(table) == "Response"), ], table$Response,
                     mtryStart=3, # Which Mtry to start with
                     ntreeTry=nTrees, # How many decision trees to build to test each Mtry
                     stepFactor=1.5, # Used to move to next Mtry to test (current + 1.5*current)
                     improve=0.0001, # Threshold of improvement if next Mtry doesn't reduce model error by this stop
                     trace=TRUE, # Print progress information
                     plot=TRUE) # Plot progress information

# Plot the output from the tuning process
plot(tuneOutput, las=1, type="o")

# Get the mtry from the tuning output
# Looks complicated but basically just taking Mtry that resulted in minimum model error
mTry <- as.integer(rownames(tuneOutput)[tuneOutput[,2] == min(tuneOutput[,2])])

# Highlight chosen Mtry on plot
points(x=mTry, y=tuneOutput[which(tuneOutput[, 1] == mTry), 2], col="red",
       pch=20)

# Run Random Forest Algorithm on FULL dataset
infoRF <- randomForest(table$Response~.,
                       data=table[, ],
                       proximity=FALSE, # Turn off similarity comparison between rows - MIGHT BE USEFUL TO CHANGE THIS?
                       mtry=mTry, # Set Mtry - taken from tuning
                       importance=TRUE, # Record the importance of the predictor variables
                       ntree=nTrees, # Number of decision trees to build
                       do.trace=FALSE, # Don't print progress information
                       keep.forest=TRUE, # Need to keep decision trees for prediction
                       norm.votes=FALSE) # Don't normalise the voting information for predictor variables

# Plot the output from the Random Forest fitting 
plot(infoRF, las=1)


# Use the model to predict the genetic distances for the test dataset
# Can we predict inter-isolate genetic distances well based only on 
# the predictor variables?
table$Predicted <- predict(infoRF, table)

# Get the Pseudo RSquared value from the fitted random forest model
# Can we explain a large amount of variation in the inter-isolate genetic distances?
# 1 is all, 0 is none
rSq <- round(infoRF$rsq[length(infoRF$rsq)], digits=2)

# Get the importance scores of each predictor variable
# Which ones are important? Have the highest influence on the model if removed
importance <- infoRF$importance



