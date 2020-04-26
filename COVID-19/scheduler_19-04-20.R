# Load the rmarkdown library
library(rmarkdown)

# Get the current date and time
timeStamp <- Sys.time()

# Note the path to webste
websiteDirectory <- file.path("~", "Desktop", "JosephCrispell.github.io", "COVID-19")

# Note the path to Rmd script
scriptsDirectory <- file.path("~", "GeneralTools", "COVID-19")

# Run the COVID-19 dahsboard Rmd file
rmarkdown::render(file.path(scriptsDirectory, "COVID-19_Scotland_01-04-20.Rmd"))

# Copy the html into my website folder
file.copy(file.path(scriptsDirectory, "COVID-19_Scotland_01-04-20.html"), websiteDirectory, overwrite=TRUE)

# Push the changes to github
setwd(file.path("~", "Desktop", "JosephCrispell.github.io"))
system("git add *")
system(paste0("git commit -m \"Daily COVID-19 dashboard update: ", timeStamp, "\""))
system("git push")
