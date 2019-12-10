# Working with genomic data

The current worksheet aims to briefly documents the steps involved in the processing and analysis of whole genome sequencing data.

## A first look at raw sequencing data

### The `FASTQ` file

Sequencing data generally comes in `fastq` format. The `fastq` file contains each of the sequencing reads with its quality information:
```
@NB501589:38:H3HTTAFXY:1:11101:16760:1040 1:N:0:GCTCATGA+NATGCAGT
GCTCGNTGGTGGTGACCGTCGCGCTGCGCGGGCTATCGCGCAGCAGCGCGATTTCGCCGACGATCATGCCCGGCAGCGCCCGAGCGATGATCGCAACACCATCGTCGCCAACATGGCTGACTTCTGCGCTACCCGACGAGATAAGCAGAAA
+
AAAAA#EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEAEEEAEAEEEAEEEEEEEAAEEEEAEEEEEEEEEEEEAEEE<AEEEA/EEEEEEAAAEE<AE
```

The above text represents the information stored for a single read produced by an illumina NextSeq machine. It represents a single fragment of an *Mycobacterium bovis* genome. Each read is represented by 4 lines of text in a `fastq` file:

Line 1 provides the following information (more info [here](https://help.basespace.illumina.com/articles/descriptive/fastq-files/)):
```
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>
```

Line 2: is a single read sequence, which corresponds to a short fragment of a genome

Line 4: contains the quality scores. Each character represents a single number based on the ASCII code of each character (more info [here](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)):

| Symbol | ASCII code | Q Score |
|--------|------------|---------|
| \!     | 33         | 0       |
| """"   | 34         | 1       |
| \#     | 35         | 2       |
| $      | 36         | 3       |
| %      | 37         | 4       |
| &      | 38         | 5       |
| '      | 39         | 6       |
| \(     | 40         | 7       |
| \)     | 41         | 8       |

### Raw read sequencing quality using fastqc

[`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is an open source piece of software that is used look at the quality of sequencing reads in a `fastq` file. To install `fastqc`, use one of the following links:

- Windows/Linux: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip
- Mac: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.dmg

You can run `fastqc` on a pair of `fastq` files either by opening them from within the fastqc program or by using the following code:
```
fastqc --threads 2 48-MBovis_R1.fastq.gz 48-MBovis_R2.fastq.gz
```
The above code will generate a selection of files. The ones we're interested in are the `.html` files. Clicking on these files will get you an output like this:

![Screenshot of output from fastqc](Images/48-MBovis_R1_fastqc.png)

The `html` files contain a series of summary plots that tell us about the quality of the sequencing data. Here is a brief description of the most useful parts of the `fastqc` output.

#### Basic statistics

The first part of the `fastqc` output provides important statistics that allow us to quickly check whether the sequencing run went ok. It is important to check the **Total Sequences**, **Sequence length**, and **%GC** values.

> QUESTION:<br>
> 1. What values would you expect for these three statistics?

For *M. bovis* genomic data, I like to have more than 250,000 read sequences, the read length is usually 150bp but that depends on the run, and the %GC is 65%.

> QUESTION:<br>
> 1. Why does the sequence length range from 35 to 151bp?

#### Per base sequence content

Each read represents a random fragment of the genomic DNA or its complement. Therefore, on average for every `A` we would expect to see a `T`, for every `C` a `G` and so on...

The **Per base sequence content** plot uses this to check the quality of the reads, where the proportion of `A`s versus `T`s or `C`s versus `G`s deviates at a position in the reads it can be considered a signal of poor quality.

Generally the accuracy of sequencing deteriorates towards the ends of the reads, so we will often see deviations at the start and end of this figure.

The **Per base sequence content** plot is useful for defining your trimming parameters.

> QUESTIONS:<br>
> 1. How much would you trim off the right?
> 2. How much off the left?

#### Per sequence GC content

Each organism's genome has a particular bias in its nucletide content. For example, *M. bovis* has a GC rich genome (lots of `C`s and `G`s). As each read represents a random fragment of the genome, and we have hundreds of thousands of them, so by counting the numbers of `C`s and `G`s in our reads we can get a good estimate of the genomes GC bias.

I find the **Per sequence GC content** plot really useful for spotting contamination. If we have accidentally sequenced a different species as well as *M. bovis*, and it has a different GC bias, this plot will look strange and our average GC content won't be 65%.

#### Further reading

Make sure and take a look through the rest of the `fastqc` output. Also, the `fastqc` help pages provides examples of [good](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) and [bad](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html) fastqc outputs that would worth taking a look at.


## Processing raw sequencing data

### Trimming sequencing reads

As we saw from the `fastqc` outputs, most sequencing data how some minor quality issues. By trimming the raw sequencing data we can improve the efficiency of our downstream analyses without impacting their accuracy.

There are loads of tools available for trimming raw sequencing data in `fastq` format. [Here](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085024) is a really nice summary of the different tools available and the benefits of trimming.

The important thing with trimming is to be conservative, if we're too harsh we might end up removing informative data. In contrast, if we're too lenient then it probably just take a wee bit longer.

I use `cutadapt`, to install it use one of the following options (more info [here](https://cutadapt.readthedocs.io/en/stable/installation.html)):

- Linux: `sudo apt install cutadapt`
- Windows: see info [here](https://www.biostars.org/p/85788/)
- Mac: `sudo easy_install cutadapt`

In the command line I can trim my raw sequence reads with the following command:

```
# Set the number of threads
NTHREADS=4

# Set your trimming parameters here
QUALITY=25 # Quality threshold used to filter poor quality reads
MINLENGTH=50 # Threshold to filter reads out that are too short
TRIMLEFTFORWARD=10 # Number of sites to remove from left of forward reads
TRIMRIGHTFORWARD=0 # Number of sites to remove from right of forward reads
TRIMLEFTREVERSE=10 # Number of sites to remove from left of forward reads
TRIMRIGHTREVERSE=0 # Number of sites to remove from right of forward reads

# Note the names of the input raw sequencing fastq files
FORWARD="48-MBovis_R1.fastq.gz"
REVERSE="48-MBovis_R2.fastq.gz"

# Name the output files
TRIMMEDFORWARD="48-MBovis_R1_trimmed.fastq.gz"
TRIMMEDREVERSE="48-MBovis_R2_trimmed.fastq.gz"

# Run cutadapt
cutadapt -b AGATCGGAAGAG -B AGATCGGAAGAG -b CTGTCTCTTATA -B CTGTCTCTTATA -o $TRIMMEDFORWARD -p $TRIMMEDREVERSE $FORWARD $REVERSE --minimum-length=$MINLENGTH --quality-cutoff=$QUALITY,$QUALITY -u $TRIMLEFTFORWARD -u -$TRIMRIGHTFORWARD -U $TRIMLEFTREVERSE -U -$TRIMRIGHTREVERSE --cores=$NTHREADS
```

The above code uses `cutadapt` to trim the forward and reverse sequence reads for a single sample (`48-Mbovis`). The above quality settings are most of the ones I commonly use. Based on the **Per base sequence content** plot from `fastqc`, I usually change the four different trimming parameters (`TRIMLEFTFORWARD`, `TRIMRIGHTFORWARD`, `TRIMLEFTREVERSE`, `TRIMRIGHTREVERSE`).

`cutadapt` is an extremely flexible, take a look at some of its additional parameters with the following command (more info [here](https://cutadapt.readthedocs.io/en/stable/guide.html)):
```
cutadapt --help
```

### Aligning the trimmed reads

We generally work with *Mycobacterium* complex species, these are haploid bacteria with conserved genomes that have little or no recombination. For *M. bovis*, these traits mean that there is little variation between samples and we can use reference genome alignment to look at the variation in our samples.

Reference genome alignment maps sequence reads against a known genome sequence (reference) that is likely to be very similar to the sample.

For *M. bovis* the [AF2122/97](https://www.ncbi.nlm.nih.gov/nuccore/LT708304) is usually used.

One of the most commonly used tools for reference genome alignment is [`bwa`](http://bio-bwa.sourceforge.net/). We can install bwa on a Linux machine using the following code:
```
sudo apt install bwa
```

Next, before we do our alignment, we'll first need to download the `fasta` file for the reference genome. You can download the `fasta` file by going to [this](https://www.ncbi.nlm.nih.gov/nuccore/LT708304) page and clicking "Send to", select "File", and then select the `FASTA` format. Put the downloaded file into a folder (I called mine `"Reference"`). We can index the reference genome `fasta` file with the following command:
```
bwa index Reference/Mbovis_LT708304-1.fasta
```

Now we can align our trimmed data for the sample (`48-MBovis`), which we produced with `cutadapt`, with the following code:
```
# Note the names of the trimmed read files
TRIMMEDFORWARD="48-MBovis_R1_trimmed.fastq.gz"
TRIMMEDREVERSE="48-MBovis_R2_trimmed.fastq.gz"

# Set the number of threads
NTHREADS=4

# Name the output file
SAMFILE="48-MBovis_aligned.sam"

# Note the path and name of reference genome
REFERENCE="Reference/Mbovis_LT708304-1.fasta"

# Run the alignment
bwa mem -t $NTHREADS $REFERENCE $TRIMMEDFORWARD $TRIMMEDREVERSE > $SAMFILE
```

The aligned data is stored in a `SAM` file (Sequence Alignment Map), which is similar in structure to the `fastq` file, except both forward and reverse reads are present in a single file and the coordinates of each read on the reference genome are recorded (more info [here](http://samtools.github.io/hts-specs/SAMv1.pdf)).

In order to view the variation that is present in our sample, we need to use a program called [samtools](http://samtools.sourceforge.net/). `samtools` is a collection tools that are used to view and process the information in a `SAM` file.

We can install samtools on a Linux machine with the following command:
```
sudo apt install samtools
```

We're going to use `samtools` to process our `SAM` file and get it ready for us to examine the variation. Run these commands:

```
# Note the name of the SAM file
SAMFILE="48-MBovis_aligned.sam"

# Set the number of threads
NTHREADS=4

# Convert SAM to BAM
BAMFILE=$NAME"_aligned.bam"
samtools view --threads $NTHREADS -b $SAMFILE > $BAMFILE

# Sort the aligned reads
SORTEDBAMFILE=$NAME"_aligned_sorted.bam"
samtools sort $BAMFILE -o $SORTEDBAMFILE --threads $NTHREADS

# Index the sorted reads
samtools index $SORTEDBAMFILE -@ $NTHREADS

# Remove duplicates
NODUPLICATES=$NAME"_aligned_sorted_rmDuplicates.bam"
samtools rmdup $SORTEDBAMFILE $NODUPLICATES
```

This is a quick test!
