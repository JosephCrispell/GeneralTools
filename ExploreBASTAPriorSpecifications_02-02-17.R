###########################
# Substitution Rate Prior #
###########################

# Open a PDF
path <- "C:/Users/Joseph Crisp/Desktop/UbuntuSharedFolder/Woodchester_CattleAndBadgers/NewAnalyses_02-06-16/BASTA/"

file <- paste(path, "PriorSpecifications_01-02-17.pdf", sep="")
pdf(file)

# Note the input for the prior ditribution
nSitesUsed <- 4049006
upper <- 0.0000003
lower <- 0.00000001
meanLog <- -15
sdLog <- 6

# Draw from the distribution
quantiles <- seq(from=lower, to=upper, by=0.000000005)
values <- dlnorm(quantiles, meanlog=meanLog, sdlog=sdLog)

values <- values / max(values)

# Convert to per genome per year rate
nSites <- 4049006
quantiles <- quantiles * nSites

# Plot the prior distribution
plot(x=quantiles, y=values, las=1, type="l",
     xlab="Substitution Rate (per genome per year)",
     main="Substitution Rate Prior Distribution",
     yaxt='n', bty='n', ylab="", ylim=c(0,1))
polygon(x=c(quantiles, quantiles[length(quantiles)], quantiles[1]),
        y=c(values, 0, 0), col=rgb(0,0,0, 0.5), border="black")
legend("topright", legend=c(paste("Lower = ", round(lower * nSites, digits=2), sep=""), 
                            paste("Upper = ", round(upper * nSites, digits=2), sep="")),
       bty="n")

################################
# Transition-Transversion Bias #
################################

# Note the input for the prior distribution
meanLog <- 0
sdLog <- 4
bound <- 10 # Ease interpretation

# Draw from the distribution
quantiles <- seq(from=0, to=bound, by=0.1)
values <- dlnorm(quantiles, meanlog=meanLog, sdlog=sdLog)
values <- values / max(values)

# Plot the prior distribution
plot(x=quantiles, y=values, las=1, type="l",
     xlab="Transition/Transversion Ratio",
     main="Transition Transversion Bias Prior Distribution",
     yaxt='n', bty='n', ylab="", ylim=c(0,1))
polygon(x=c(quantiles, quantiles[length(quantiles)], quantiles[1]),
        y=c(values, 0, 0), col=rgb(0,0,0, 0.5), border="black")
legend("topright", legend=c("Lower = 0", 
                            expression(paste("Upper = ", infinity, sep=""))),
       bty="n")

#########################
# State Transition Rate #
#########################

# Note the input for the prior distribution
mean <- 0.01
bound <- 0.05 # Ease interpretation

# Draw from the distribution
quantiles <- seq(from=0, to=bound, by=0.001)
values <- dexp(quantiles, rate=1/mean)
values <- values / max(values)

# Plot the prior distribution
plot(x=quantiles, y=values, las=1, type="l",
     xlab="Rate (Events/Year)",
     main="State Transition Rate Prior Distribution",
     yaxt='n', bty='n', ylab="", ylim=c(0,1))
polygon(x=c(quantiles, quantiles[length(quantiles)], quantiles[1]),
        y=c(values, 0, 0), col=rgb(0,0,0, 0.5), border="black")
legend("topright", legend=c("Lower = 0", 
                            expression(paste("Upper = ", infinity, sep=""))),
       bty="n")

####################
# Population Sizes #
####################

# Note the input for the prior distribution
meanLog <- 0
sdLog <- 1
bound <- 50 # Ease interpretation

# Draw from the distribution
quantiles <- seq(from=0, to=bound, by=0.1)
values <- dlnorm(quantiles, meanlog=meanLog, sdlog=sdLog)
values <- values / max(values)

# Plot the prior distribution
plot(x=quantiles, y=values, las=1, type="l",
     xlab="Effective Size",
     main="Population Size Prior Distribution",
     yaxt='n', bty='n', ylab="", ylim=c(0,1))
polygon(x=c(quantiles, quantiles[length(quantiles)], quantiles[1]),
        y=c(values, 0, 0), col=rgb(0,0,0, 0.5), border="black")
legend("topright", legend=c("Lower = 0", 
                            expression(paste("Upper = ", infinity, sep=""))),
       bty="n")

dev.off()
