#### Preparation ####

# Load libraries
library(xlsx)

# Set the path
path <- file.path("~", "Desktop")

# Get today's date
date <- format(Sys.Date(), "%d-%m-%y")

#### Read in the data ####

# Read in the Author ID table
authorIDs <- read.xlsx(file.path(path, "637075988602021131_pubmedresult-14.xml.xlsx"), sheetName = "Author")

# Read in the Affiliation info
affiliations <- read.xlsx(file.path(path, "637075988602021131_pubmedresult-14.xml.xlsx"), sheetName = "AffiliationInfo")

#### Identify UCD authors ####

# Find author IDs for those affiliated with UCD
ucdAuthorIDs <- affiliations[grepl(affiliations$Affiliation, pattern="University College Dublin"), "Author_Id"]

# Get first and surnames of authors affiliated with UCD
authorNames <- authorIDs[authorIDs$Author_Id %in% ucdAuthorIDs, c("ForeName", "LastName")]
nrow(authorNames)
