#### Preparation ####

# Load the required libraries
library(ape) # Reading in FASTA
library(lubridate) # Converting to decimal dates

# Set the path
path <- file.path("~", "Desktop", "BuildingBASTAXML", "BASTA_equal_relaxed_25-02-20")

# Set the date
date <- format(Sys.Date(), "%d-%m-%y")