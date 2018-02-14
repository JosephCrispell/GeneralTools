#### Packages ####

library(randomForest)

#### Read in data ####

table <- read.table("fileName", header=TRUE, stringsAsFactors=FALSE, sep=",")

#### Random Forest model Fitting ####

# Split the data into train and test sets
trainProp <- 0.5
trainRows <- sample(x=1:nrow(table),
                    size=floor(trainProp * nrow(table)), replace=FALSE)

# Define the columns of the table we're interested in
cols <- 2:26
nTrees <- 1000

# Run the Random Forest tuning algorithm to find the optimal mtry value
tuneOutput <- runRandomForestTuning(table[trainRows, ], 3, nTrees, cols)
plot(tuneOutput, las=1, type="o")

# Calculate the mtry from the tuning output
mTry <- as.integer(rownames(tuneOutput)[tuneOutput[,2] == min(tuneOutput[,2])])
points(x=mTry, y=tuneOutput[which(tuneOutput[, 1] == mTry), 2], col="red",
       pch=20)

# Run Random Forest Algorithm
nTrees <- 1000
infoRF <- randomForest(table[trainRows, "GeneticDistance"]~., data=table[trainRows, cols],
                       proximity=FALSE, mtry=mTry, importance=TRUE,
                       ntree=nTrees, do.trace=FALSE, keep.forest=TRUE,
                       norm.votes=FALSE)
plot(infoRF, las=1)

#### Plot the results ####

# Add the predicted values onto the table
predictions <- predict(infoRF, table[-trainRows, cols])

# Get the Pseudo RSquared value
rSq <- round(infoRF$rsq[length(infoRF$rsq)], digits=2)

# Predicted Versus Actual
smoothScatter(predictions, table[-trainRows, "GeneticDistance"],
              main="Actual vs. Predicted",
              cex.main=1,
              xlab="Predicted",
              ylab="Actual",
              nrpoints=0, las=1)
abline(lm(table$GeneticDistance ~ table$predictions), col="red")
correlation <- round(cor(table[-trainRows, "GeneticDistance"], predictions), digits=2)
legend("topleft", legend=c(paste("corr =", correlation), paste("Rsq =", rSq)), bty="n", cex = 1)

# Variable importance
plotVariableImportance(infoRF=infoRF, colToUse="%IncMSE")



#### FUNCTIONS ####

plotVariableImportance <- function(infoRF, colToUse){
  
  # Get the variable importance from the RF model
  variableImportance <- as.data.frame(infoRF$importance)
  
  # Order the table by importance
  variableImportance <- variableImportance[order(variableImportance[, colToUse],
                                                 decreasing=FALSE), ]
  
  # Transpose the table
  transpose <- as.matrix(t(variableImportance))
  

  plot <- barplot(transpose[-2,], horiz=TRUE, beside=TRUE,
                  main="Variable Importance",
                  col.axis="white")
}


runRandomForestTuning <- function(inputTable, initialMtry, nTrees, cols){
  
  tuneOutput <- tuneRF(inputTable[, cols], inputTable$Genetic, mtryStart=initialMtry,
                       ntreeTry=nTrees, stepFactor=1.5, improve=0.0001, trace=TRUE,
                       plot=TRUE)
  
  return(tuneOutput)
}

