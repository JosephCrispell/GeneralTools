# Load libraries
library(reticulate) # Running python code

# Set to use python3
use_python("/usr/local/bin/python3", required=TRUE)
py_config() # To check

# Load the googletrans python package
py_run_string("from googletrans import Translator")

# Load the python translation code
py_run_string("translator = Translator()")

# Source the python translation code
source_python(file.path("", "Users", "josephcrispell", "Desktop", "GeneralTools", "TranslateUsingGoogleTrans.py"))

# Translate some text
text <- "Traduire du texte à l'aide de googletrans à partir de R."
translation <- translate(text, 'en')

# Create typing animation
translationAnimation(text, translation)

# FUNCTIONS
translationAnimation <- function(french, english){
  
  # Type out the french
  typeText(french)
  
  # Delete the french
  typeText(french, backward=TRUE)
  
  # Type out the english
  typeText(english)
}

typeText <- function(text, backward=FALSE, wait=0.1, nCharToPad=300){
  
  # Define the pad to use before and after text
  pad=paste(rep(" ", 100), collapse="")
  
  # Split the text into characters
  characters <- strsplit(text, split="")[[1]]
  
  # Print text to screen if going backwards
  if(backward){
    cat(paste0("\r", pad, text, pad))
  }
  
  # Print each character one at a time
  for(i in seq_along(characters)){
    
    # Wait set period
    Sys.sleep(wait)
    
    # Not the characters to print
    charactersToPrint <- characters[1:i]
    if(backward){
      charactersToPrint <- characters[1:(nchar(text)-i)]
    }
    charactersToPrint <- c(charactersToPrint, rep(" ", nchar(text)-length(charactersToPrint)))
    
    # Print next character
    cat(paste0("\r", pad, paste(charactersToPrint, collapse=""), pad))
  }
  
  # Clear line if printing backwards
  if(backward){
    cat(paste0("\r", pad, paste(rep(" ", nchar(text)), collapse=""), pad))
  }
}
