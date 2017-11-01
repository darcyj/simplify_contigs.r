# simplify_contigs.r
Simplifies output from genome/metagenome assembly by combining contigs with redundant regions
by Jack Darcy, Oct 2017

My metagenome assembly pipeline is as follows:
1. Aggressively quality filter paired-end reads using vsearch (https://github.com/torognes/vsearch)
2. Assemble filtered reads together using SPAdes (http://bioinf.spbau.ru/en/spades)
3. Use MaxBin (https://sourceforge.net/projects/maxbin/) to cluster contigs together into bins, which hopefully contain contigs from the same original genome (or at least enrich for contigs from the same genome!). 

However, I have found that many contigs in the same bin share large overlapping regions, and likely should be merged into the same contig. While this kind of merging may not be desirable for analyses of syntony, if one is only interested in the functional capacity of a bin then it is desirable to eliminate redundancy from the bin before downstream analysis (annotation, etc), and this process has the benefit of producing longer contigs, which is great. The program CAP3 *kind of* does this (http://seq.cs.iastate.edu/cap3.html), but the output is cryptic, there are too many options, and it's old as hell. It was written for Sanger data. I want to be sure that I'm only combining sequences if they are 100% identical in the overlapping region. Illumina data are AWESOME in terms of quality, even moreso if you aggressively filter. Plus, the added confidence you get when a contig has a depth (number of stacked-up sequences) above 3 or 4. That's what this script relies on. 

Usage:

simplify_contigs.r -i contigs.fasta -o simplified_contigs.fasta -m 20 -w 50

-i/--input: your un-simplified, freshly-binned or freshly-assembled contigs in SEQUENTIAL FASTA FORMAT. NOT INTERLEAVED.

-o/--output: filepath for simplified contigs

-m/--maxiters: number of times this script should iteratively combine sequences. Rarely goes over 4. Default=20

-w/--wordlength: number of characters to use for sequence-sequence matching. this is also the minimum overlap to merge two sequences. Default=50

-c/--skip_nested_check: Usage depends on the assembler. I know SPAdes will never output two contigs where one is perfectly nested within the other, so the default is not to check for that condition before combining. But other assemblers may be different. 
