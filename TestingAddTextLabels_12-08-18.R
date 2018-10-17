#library(addTextLabels)

# Create some random points
n <- 45
testLabels <- c("short", "mediummm", "looooonnnnnnngggggg", "0090292002", "9", "A Different label")
coords <- data.frame(X=runif(n), Y=runif(n, min=0, max=100), Name=sample(testLabels, size=n, replace=TRUE),
                     stringsAsFactors = FALSE)

# Plot them without labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")

pdf("/home/josephcrispell/Desktop/test.pdf")

# With potentially overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
text(coords$X, coords$Y, labels=coords$Name, xpd=TRUE)

# Plot them with non-overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
addTextLabels(coords$X, coords$Y, coords$Name, cex=0.9, col.label="black")

# Plot them with non-overlapping labels
plot(x=coords$X, y=coords$Y, pch=19, bty="n", xaxt="n", yaxt="n", col="red", xlab="X", ylab="Y")
addTextLabels(coords$X, coords$Y, coords$Name, cex=0.9, col.background=rgb(0,0,0, 0.75), col.label="white")

dev.off()


x = c(0.8846, 1.1554, 0.9317, 0.9703, 0.9053, 0.9454, 1.0146, 0.9012, 
      0.9055, 1.3307)
y = c(0.9828, 1.0329, 0.931, 1.3794, 0.9273, 0.9605, 1.0259, 0.9542, 
      0.9717, 0.9357)
ShortSci = c("MotAlb", "PruMod", "EriRub", "LusMeg", "PhoOch", "PhoPho", 
             "SaxRub", "TurMer", "TurPil", "TurPhi")

# basic plot
plot(x, y, asp=1)
abline(h = 1, col = "green")
abline(v = 1, col = "green")

addTextLabels(x, y, ShortSci, cex=0.9, col.background=rgb(0,0,0, 0.75), col.label="white")
addTextLabels(x, y, ShortSci, cex=0.9, col.background=rgb(0,0,0, 0.75), col.label="white", keepLabelsInside = FALSE)
