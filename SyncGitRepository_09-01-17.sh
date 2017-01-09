#!/bin/sh
#!/bin/bash

# Sync local UNIX repository with remote git server
# Author: Joseph Crisp

# RUN SCRIPT USING bash script.sh (in Directory containing Raw Reads)
# Note create permissions first ---> chmod +x script.sh

# Command Line Structure:
# bash SyncGitRepository.sh username password

###################
# Necessary Input #
###################

####Github account details - from command line
USERNAME=$1
PASSWORD=$2

#####################
# Update Repository #
#####################

#### Find untracked (new) files
NEWFILES=`git ls-files --others --exclude-standard`

#### Add new files to git repository
for FILE in ${NEWFILES[@]}; 
do
    git add $FILE
done