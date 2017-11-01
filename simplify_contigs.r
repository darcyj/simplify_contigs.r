#!/usr/bin/env Rscript

## set up environment
	source("/ResearchSoftware/jld_scripts/simplify_contigs_functions.r")
	suppressPackageStartupMessages(require(optparse))

## make options for script
	option_list <- list(
		make_option(c("-i", "--input"), action="store", default=NA, type='character',
			help="Input contigs filepath, in sequential (NOT INTERLEAVED) fasta format."), 
		make_option(c("-o", "--output"), action="store", default=NA, type='character',
			help="Output filepath for simplified contigs in fasta format."),		
		make_option(c("-m", "--maxiters"), action="store", default=20, type='integer',
			help="Maximum number of iterations for contig joining."),
		make_option(c("-w", "--wordlength"), action="store", default=50, type='integer',
			help="Length of DNA sequence to use for match searching."),
		make_option(c("-c", "--skip_nested_check"), action="store_true", default=FALSE, type='logical',
			help="Skips the check for nested contigs in iteration 1.")
	)
	# parse option arguments
	opt <- parse_args(OptionParser(option_list=option_list))
	input_fp <- opt$input
	output_fp <- opt$output
	word_length <- opt$wordlength
	iters <- opt$maxiters
	skip_nested_check <- opt$skip_nested_check
	
## read in fasta file as a list of strings
	print("Reading in and formatting data:")
	seqs <- scan(input_fp, what="", sep="\n")
	seq_names <- removefastacarat(simplify2array(seqs[1:length(seqs) %% 2 == 1]))
	seqs <- simplify2array(seqs[1:length(seqs) %% 2 == 0])
	seq_lengths <- simplify2array(lapply(X=seqs, FUN=nchar))


	# convert to data frame for easy management
	seqs_df <- data.frame(name=seq_names, length=seq_lengths, seq=seqs, stringsAsFactors=F)
	rm(seqs, seq_names, seq_lengths)
	print(paste("...done. Read", nrow(seqs_df), "sequences."))

## store number of input sequences
	n_seqs_input <- nrow(seqs_df)
	# initialize n_seqs_consumed to an absurd number
	# it will be overwritten within the while loop, but
	# it needs to start high so that the loop actually runs
	n_seqs_consumed <- 8008135

## iteratively simplify contigs
	i <- 0
	while(n_seqs_consumed > 0){
	
		i <- i + 1

		## print iteration message:
		print(paste("Beginning iteration", i))
		
		## sort by length
		print("...Sorting sequences:")
		seqs_df <- seqs_df[order(seqs_df$length, decreasing=T), ]
		print("......done.")

		
		## store number of sequences before iteration
		n_seqs_before <- nrow(seqs_df)
		
		## find matches and combine 'em
		print("...Finding and combining matches:")
		seqs_df <- simplify_contigs(seqs_df)
		print("......done.")
		
		## calculate number of sequences consumed
		n_seqs_consumed <- n_seqs_before - nrow(seqs_df)

		## break iterative loop if no sequences consumed.
		if(n_seqs_consumed == 0){
			print("...This iteration could not join any contigs together. Terminating.")
			break()
		}else{
			print(paste("...This iteration consumed", n_seqs_consumed, "contigs."))
		}
		
		## break iterative loop if maxiters has been reached
		if(i == iters){
			print("Maximum specified iterations reached. Terminating.")
			print("This means the program did NOT run to completion.")
			break()
		}
		
	}	
	
## After assembly, check for nested sequences
	## this requires sorting: final sort
	print("Sorting sequences final time:")
	seqs_df <- seqs_df[order(seqs_df$length, decreasing=T), ]
	print("...done.")

	## for each contig, check if it is 100% contained within another. if it is, get rid of it
	if(skip_nested_check == F){
		print("Removing nested contigs:")
		seqs_nested <- find_nested_seqs(seqs_df)
		seqs_df <- seqs_df[ !seqs_nested, ]
		print(paste("...done. Removed", sum(seqs_nested), "redundant contigs."))
	}

## State progress made, and write out fasta file
	print(paste(n_seqs_input - nrow(seqs_df), "contigs consumed by iterative assembly."))
	print(paste("Writing output fasta file to", output_fp, ":"))
	write_seqs_df_to_fasta(seqs_df, output_fp)







