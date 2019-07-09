#!/bin/sh
#!/bin/bash

# Convert all PDFs in current directory to PNG files
# Author: Joseph Crispell

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

# Command Line Structure:
# bash ConvertPDFToPNG_DATE.sh commitMessage

#############################
# Convert PDFs to PNG files #
#############################

# Get a list of all pdf files in current directory
PDFS=(`ls | grep ".pdf$"`)

# Check if PDF file provided
if [ "$#" == 1 ]
then
	PDFS=($1)
fi

# Examine each PDF
for (( PDFINDEX=0; PDFINDEX<${#PDFS[@]}; PDFINDEX+=1 ))
do
#
	# Get the current file
	FILE=${PDFS[$PDFINDEX]}

	# Count number of pages in current pdf
	NPAGES=`pdfinfo $FILE | grep "Pages" | awk '{ split($0, array, " +"); print array[2] }'`

	# Print progress
	echo -e "\e[0;34m Beginning to convert PDF: \e[0m"$FILE"\e[0;34m N. pages: \e[0m"$NPAGES

	# Examine every page of current PDF
	for (( PAGEINDEX=0; PAGEINDEX<$NPAGES; PAGEINDEX+=1 ))
	do

		# Build name for PNG file
		PNG=${FILE::-4}"_"$PAGEINDEX".png"
		
		# Convert the current page of the current PDF to PNG
		convert -density 300 $FILE[$PAGEINDEX] -quality 100 $PNG

		# Print progress
		echo -e "	Created PNG: "$PNG
	done
done

# Print progress
echo -e "\e[0;32m Finished. \e[0m"
