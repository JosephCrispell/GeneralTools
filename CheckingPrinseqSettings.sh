#!/bin/sh
#!/bin/bash

# Changing settings for Prinseq

# Command Line Structure:
# bash CheckingPrinseqSettings.sh R1.fastq R2.fastq

# Trimming - Needs: ~/prinseq-lite-0.20.4/prinseq-lite.pl installed in this location
LENGTH=50			# The minimum length of read to be accepted
MEANQUAL=20			# Filter sequence if mean quality score below x
TRIML=20			# Trim sequence at the 5' end by x positions
TRIMR=5			# Trim sequence at the 3' end by x positions
TRIMQUALL=20		# Trim sequence by quality score from the 5' end with this threshold
TRIMQUALR=20		# Trim sequence by quality score from the 3' end with this threshold score
TRIMTYPE="mean"		# Type of quality score calculation to use [min, mean, max, sum]
WINDSIZE=10			# The sliding window size used to calculate quality score by type
TRIMLTAIL=5			# Trim poly A/T > X length at 5' end
TRIMRTAIL=5			# Trim poly A/T > X length at 3' end

echo -e "\e[0;34m Beginning Read Trimming... \e[0m"

perl ~/prinseq-lite-0.20.4/prinseq-lite.pl -fastq $1 -fastq2 $2 -min_len $LENGTH -min_qual_mean $MEANQUAL -trim_left $TRIML -trim_right $TRIMR -trim_qual_left $TRIMQUALL -trim_qual_right $TRIMQUALR -trim_qual_type $TRIMTYPE -trim_qual_window $WINDSIZE  

echo -e "\e[0;31m Finished. \e[0m"