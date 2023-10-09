# Population Genomics lab notebook  

## Author: Maryam Nouri-Aiin
### Affiliation:  PSS
### E-mail contact: mnouriai@uvm.edu


### Start Date: 09-11-2023
### End Date: TBD
### Project Descriptions:   Notebook to document the bioinformatic of Population Genomics


```{css my-header-colors, echo = FALSE}
.page-header {
    background-image: linear-gradient(50deg,blue, orange);
}
```



# Table of Contents:   

* [Entry 1: 2023-09-11](#id-section1)
* [Entry 2: 2023-09-13](#id-section2)
* [Entry 3: 2023-XX-XX](#id-section3)


------    
<div id='id-section1'/>   

### Entry 1: 2023-09-11.   
- Red spruce study system and exome capture file
- Structure of fastq file (DNA sequnces, plus the Qscores)
- Using program FastQC, analyzed the sequncing run for one file
------    
<div id='id-section2'/>   


### Entry 2: 2023-09-13.  
- FastQC results discussion, good qulaity sequence data for most reads
- Initial 5 bp or so had more variable base frequncies, The very end of reads had slightly lower @-score
- Based on this we set up an analysis to trim the reads, using 'fastq' program
- bash script run  'fast.sh'
- looked at html files produced by 'fastq' and compared pre- and post-trimming -- things looked good!
- We ended the day setting up our read mapping of the trimmed and cleaned read using 'tmux'


------    
<div id='id-section3'/>   


### Entry 3: 2023-XX-XX.




<div style="background-color:#F8E0E0; padding:10px;">

### **Details 
2. Here’s our “pipeline”
Visualize the quality of raw data (Program: FastQC)

Clean raw data (Program: Trimmomatic)

Visualize the quality of cleaned data (Program: FastQC)

Calculate #’s of cleaned, high quality reads going into mapping

We’ll then use these cleaned reads to align to the reference genome next time so that we can start estimating genomic diversity and population structure.

3.-5. Visualize, Clean, and Visualize again
Whenever you get a new batch of NGS data, the first step is to look at the data quality of coming off the sequencer and see if we notice any problems with base quality, sequence length, PCR duplicates, or adapter contamination. A lot of this info is stored in the raw data files you get from the core lab after sequencing, which are in “fastq” format.

The fastq files for our project are stored in this path: /netfiles/ecogen/PopulationGenomics/fastq/red_spruce

cd over there and ll to see the files. There should be 190 fastq files – 2 for each of the 95 samples (2 files/sample because these are paired-end reads, and each sample gets a file with the forward reads (R1) and another with the reverse reads (R2)).

The naming convention for our data is: <PopCode>_<RowID>_<ColumnID>_<ReadDirection>.fast.gz

Together, <PopCode>_<RowID>_<ColumnID> define the unique individual ID for each DNA sample, and there should be 2 files per sample (and R1 and an R2)

So…what is a .fastq file anyway?
A fastq file is the standard sequence data format for NGS. It contains the sequence of the read itself, the corresponding quality scores for each base, and some meta-data about the read.

The files are big (typically many Gb compressed), so we can’t open them completely. Instead, we can peek inside the file using head. But size these files are compressed (note the .gz ending in the filenames), and we want them to stay compressed while we peek. Bash has a solution to that called zcat. This lets us look at the .gz file without uncompressing it all the way. Let’s peek inside a file:

zcat 2505_9_C_R2.fastq.gz | head -n 4

@A00354:455:HYG3FDSXY:1:1101:3893:1031 2:N:0:CATCAAGT+TACTCCTT
GTGGAAAATCAAAACCCTAATGCTGAAAGGAATCCAAATCAAATAAATATTTTCACCGACCTGTTTCGATGCCAGAATTGTCTGCGCAGAAGACTCGTGAAATTTCGCCAGCAGGTAAAATTAAAAGGCTAGAATTAACCGCTGAAATGGA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:F:F
Note: zcat lets us open a .gz (gzipped) file; we then “pipe” | this output from zcat to the head command and print just the top 4 lines -n4

The fastq file format has 4 lines for each read:

Line	Description
1	Always begins with ‘@’ and then information about the read
2	The actual DNA sequence
3	Always begins with a ‘+’ and sometimes the same info in line 1
4	A string of characters which represent the quality scores; always has same number of characters as line 2
Here’s a useful reference for understanding Quality (Phred) scores. If P is the probability that a base call is an error, then:

Q = -10*log10(P)

So:

Phred Quality Score	Probability of incorrect base call	Base call accuracy
10	1 in 10	90%
20	1 in 100	99%
30	1 in 1000	99.9%
40	1 in 10,000	99.99%
The Phred Q score is translated to ASCII characters so that a two digit number can be represented by a single character.

 Quality encoding: !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHI
                   |         |         |         |         |
    Quality score: 0........10........20........30........40   
What kind of characters do you want to see in your quality score?

Visualize using FastQC
We’re going to use the program FastQC (already installed on our server). FastQC looks at the quality collectively across all reads in a sample.

First, let’s cd back to our home directories ~/ and set up some new folders to store our work. We’ll make 3 directories to store our data, scripts, and results:

mkdir mydata/
mkdir myscripts/
mkdir myresults/
Then let’s cd into the myresults/ folder then use pwd to prove to yourself that you’re in the myresults/ folder within your home directory. It should look like this (but with your home directory info instead of mine):

[kellrlab@ecogen myresults]$ pwd
/users/k/e/kellrlab/myresults
Now within myresults/ let’s make another folder called fastqc/ to hold the outputs from our QC analysis. Do that on your own, just like we did above, then cd into the fastqc/ folder and type pwd again to prove to yourself you did it right.

Now, we’re ready to run FastQC to look at the quality of our sequencing. The basic command is like this:

fastqc filename.fastq.gz -o outputdirectory/
This will generate an .html output file for each input file you’ve run.

Once you’ve got results, let’s use Filezilla to transfer the folder ~/myresults/fastqc/ over to your laptop and into the results folder in your github repo. Once the file transfer is complete, go to where you saved your files on your laptop and try double-clicking one of the html outputs. It should open with a web browser.

How does the quality look?

Since we made some changes (added files) to our github repo, we should practice committing these and then pushing to GitHub!
</div>


