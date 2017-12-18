#!/bin/sh
#!/bin/bash

# Sync local UNIX repository with remote git server
# Author: Joseph Crispell

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

# Command Line Structure:
# bash SyncGitRepository.sh commitMessage

# Note need to assign your email and username to github
# git config --global user.email "you@example.com"
# git config --global user.name "Your Name"

######################
# Set commit message #
######################

COMMIT="GeneralUpdate"
if [ ! $1 == "" ]
then
	COMMIT=$1
fi

#####################
# Update Repository #
#####################

#### Find untracked (new) files
NEWFILES=`git ls-files --others --exclude-standard`

#### Add new files to git repository
for FILE in ${NEWFILES[@]}; 
do
    git add $FILE
	echo "Added "$FILE" to git repository."
done

#### Remove any files that have been deleted
git add . -A

#### Commit changes to repository
TODAY=`date +%d-%m-%Y`
git commit -m $COMMIT" "$TODAY

#### Push any changes in repository to git server
git push
