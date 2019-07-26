#!/bin/bash

cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o trimmed_forward.fastq.gz -p trimmed_reverse.fastq.gz input_forward.fq.gz input_reverse.fq.gz

#-b adapter can be ligated to 3' or 5' end on forward read
#-B as above but for reverse

#-u/U LENGTH --cut=LENGTH Remove X bases from beginning (if +ve) or end (if -ve) of reads
#--trim-n Trim Ns on ends of reads
#--minimum-length Discard short reads
#--max-n Discard reads with more than X bases. If number between 0 and 1 then intepreted as proportion
#-q [5'cutoff,3'cutoff] --quality-cutoff=[5'cutoff,3'cutoff] Trim low-quality bases from 5' and/or 3' ends of each read

#The -A/-G/-B/-U options work like their -a/-b/-g/-u counterparts, but
#are applied to the second read in each pair.
#-A ADAPTER          3' adapter to be removed from second read in a pair.
#-G ADAPTER          5' adapter to be removed from second read in a pair.
#-B ADAPTER          5'/3 adapter to be removed from second read in a pair.
#-U LENGTH           Remove LENGTH bases from second read in a pair (see --cut).

#--cores=10

cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o trimmed_forward.fastq.gz -p trimmed_reverse.fastq.gz input_forward.fq.gz input_reverse.fq.gz --minimum-length=50 --trim-n --max-n=0.25 --quality-cutoff=25,25 -u 10 -u -10 -U 10 -U -10 --cores=10

	# Define Prinseq settings - USE THIS FORMAT FOR INPUT FILE IF USING ONE
	LENGTH=50			# The minimum length of read to be accepted
	MEANQUAL=20			# Filter sequence if mean quality score below x
	TRIML=20			# Trim sequence at the 5' end by x positions
	TRIMR=5				# Trim sequence at the 3' end by x positions
	TRIMQUALL=20		# Trim sequence by quality score from the 5' end with this threshold
	TRIMQUALR=20		# Trim sequence by quality score from the 3' end with this threshold score
	TRIMTYPE="mean"		# Type of quality score calculation to use [min, mean, max, sum]
	WINDSIZE=10			# The sliding window size used to calculate quality score by type
	TRIMLTAIL=5			# Trim poly A/T > X length at 5' end
	TRIMRTAIL=5			# Trim poly A/T > X length at 3' end

	LENGTH=50           # The minimum length of read to be accepted
	QUAL=25             # Trim low-quality bases from 5' and 3' ends
	PROPN=0.5           # Discard reads with more than this proportion of Ns
	TRIML=20            # Trim sequence at the 5' end by x positions
	TRIMR=5             # Trim sequence at the 3' end by x positions
