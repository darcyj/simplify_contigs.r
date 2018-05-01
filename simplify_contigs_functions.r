# remove > from start of sequence names
removefastacarat <- function(seqids){
	seqids <- substr(seqids, start=2, stop=nchar(seqids))
	return(seqids)
}


# functions to reverse-compliment a sequence
complimnet <- function(x){
	out <- "-"
	if(x == "C"){out <- "G"}
	if(x == "G"){out <- "C"}
	if(x == "A"){out <- "T"}
	if(x == "T"){out <- "A"}
	return(out)
}
rc <- function(x){
	if(length(x) > 0){
		for(i in 1:length(x)){
			tmp_i <- rev(unlist(strsplit(x[i], split="")))
			x[i] <- paste(sapply(tmp_i, FUN=complimnet), collapse="")
		}
	}
	return(x)
}

# function to look for a sequence inside another sequence,
	# even when the search term is > 2559 characters
strXinstrY <- function(y, x){
	answer <- FALSE
	if (nchar(x) < 2559){
		answer <- grepl(x, y)
	}else{
		search_word <- substr(x, 0, 2559)
		start_positions <- unlist(gregexpr(search_word, y, fixed=T))
		if(sum(start_positions) > -1){
			for(i in 1:length(start_positions)){
				sub_start_i <- start_positions[i]
				sub_stop_i <- start_positions[i] + nchar(x) - 1
				if(sub_stop_i <= nchar(y)){
					y_i <- substr(y, start=sub_start_i, stop=sub_stop_i)
					answer_i <- x == y_i
					if(answer_i == T){
						answer <- TRUE
						break()					
					}
				}
			}
		}
	}
	return(answer)
	# test code:
	# x <- paste(rep("abcdefghijklmnopqrst", 500), collapse="")
	# nchar(x) 
	# y <- paste(	paste(rep("aldiahcnmaldoalskdfh", 500), collapse=""), 	x, 	paste(rep("aklsjdhflaiwuefalksdjhfais", 100), collapse=""), 	collapse="")
	# strXinstrY(x, y)
}

grepl_long <- function(pattern, x){
	return(as.vector(sapply(X=x, FUN=strXinstrY, x=pattern)))
}

find_matches <- function(in_seq, seqs2search, exclude=0, search_word){
	# find FWD RS matches
	match_positions <- regexpr(search_word, text=seqs2search, fixed=T)
	match_positions[exclude] <- -1
	
	# find RC RS matches
	search_word_rc <- rc(search_word)
	match_positions_rc <- regexpr(search_word_rc, text=seqs2search, fixed=T)
	match_positions_rc[exclude] <- -1
	
	# count the number of matching sequences found
	n_matches <- sum(match_positions_rc > -1) + sum(match_positions > -1)
	
	output <- rep("miss", length(seqs2search))
	output[match_positions > -1] <- "FWD_hit"
	output[match_positions_rc > -1] <- "RC_hit"

	return(output)
}

testalign <- function(seq1, seq2, search_word, return_s1_leftside=F){
	# word positions in seq1
	wordstarts_1 <- unlist(gregexpr(search_word, seq1, fixed=T))
	
	# word positions in seq2
	wordstarts_2 <- unlist(gregexpr(search_word, seq2, fixed=T))
	
	# table of combinations of wordstarts
	wordstarts_table <- expand.grid(seq1=wordstarts_1, seq2=wordstarts_2)
	
	# calculate overlap length for each row in table
	seqlengths <- c(nchar(seq1), nchar(seq2))
	overlap <- rep(0, nrow(wordstarts_table))
	for(i in 1:nrow(wordstarts_table)){
		left_nonoverlaps <- wordstarts_table[i,] - 1
		overlap[i] <- min(seqlengths - left_nonoverlaps)
	}
	wordstarts_table <- data.frame(wordstarts_table, overlap)
	
	# sort by overlaps, check longest first
	wordstarts_table <- wordstarts_table[order(wordstarts_table$overlap, decreasing=T), ]

	# check each combo
	output <- FALSE
	for(i in 1:nrow(wordstarts_table)){
		left_nonoverlaps <- wordstarts_table[i,] - 1
		seq1_overlap <- substr(seq1, start=wordstarts_table[[1]][i] - min(left_nonoverlaps), stop=(wordstarts_table[[1]][i] + wordstarts_table$overlap[i] - 1))
		seq2_overlap <- substr(seq2, start=wordstarts_table[[2]][i] - min(left_nonoverlaps), stop=(wordstarts_table[[2]][i] + wordstarts_table$overlap[i] - 1))
		if (seq1_overlap == seq2_overlap){
			if(return_s1_leftside){
				output <- substr(seq1, start=1, stop= (wordstarts_table[[1]][i] - 1))
			}else if (sum(left_nonoverlaps) > 0){ # even if it's a good alignment, it's bad if there's no overhang.
				output <- TRUE
			}
			break()		
		}else if(return_s1_leftside){
			output <- "FAIL"		
		}
	}
	
	return(output)

}

unique_combos <- function(x, y, preservecolnames=F){
	df1 <- expand.grid(x, y)
	df1 <- df1[!duplicated(t(apply(df1, 1, sort))),]
	if(preservecolnames){colnames(df1) <- c(deparse(substitute(x)), deparse(substitute(y)))}
	return(df1)
}


find_nested_seqs <- function(seqs_df){
	# create output, initialize to FALSEs
	seq_within <- rep(F, nrow(seqs_df))
	
	# start up progress bar
	pb <- txtProgressBar(min = 0, max = nrow(seqs_df), style = 3)
	
	# find nested sequences
	for(i in 2:nrow(seqs_df)){
		seq_i <- seqs_df$seq[i]
		larger_tf <- seqs_df$length >= nchar(seq_i)
		larger_tf[i] <- F
		n_matches <- sum(grepl_long(pattern=seq_i, x=seqs_df$seq[larger_tf]))
		
		# if no matches in FWD, check REV:
		if (n_matches == 0){
			seq_i <- rc(seq_i)
			n_matches <- sum(grepl_long(pattern=seq_i, x=seqs_df$seq[larger_tf]))
		}
		
		if(n_matches > 0){seq_within[i] <- TRUE}
		
		# update progress bar
		setTxtProgressBar(pb, i)
		
	}
	return(seq_within)
}


extend_contig <- function(seqs_df, contig_index, word_length=50, direction="left"){
	
	i <- contig_index
	seq_i <- seqs_df$seq[i]
	
	# go right by RC-ing seq_i
	if(direction == "right"){seq_i <- rc(seq_i)}

	# left side matching - only search against sequences SHORTER than seq_i
	# this is done so that when matches are removed, the script can't accidentally skip sequences.
	search_word <- substr(seq_i, start=1, stop=word_length)
	matches_i <- find_matches(in_seq=seq_i, seqs2search=seqs_df$seq, exclude=1:i, search_word)
	
	#sum( matches_i != "miss")
		
	if( ! all(matches_i == "miss") ){
		match_seqs <- seqs_df$seq[matches_i == "FWD_hit"]
		match_seqs_indices <- which(matches_i == "FWD_hit")
		
		# rc the rc matches
		match_seqs <- c(match_seqs, rc(seqs_df$seq[matches_i == "RC_hit"]))
		match_seqs_indices <- c(match_seqs_indices, which(matches_i == "RC_hit"))

		
		# order matches by length, so longest is first
		match_seqs_order <- order(nchar(match_seqs), decreasing=T)
		match_seqs<- match_seqs[match_seqs_order]
		match_seqs_indices <- match_seqs_indices[match_seqs_order]
		
		# for each match, see if it's a good match
		match_good <- rep(F, length(match_seqs_indices))
		for(j in 1:length(match_seqs_indices)){
			match_good[j] <- testalign(match_seqs[j], seq_i, search_word, return_s1_leftside=FALSE)
		}
		
		# drop bad matches
		match_seqs <- match_seqs[match_good]
		match_seqs_indices <- match_seqs_indices[match_good]
		
		# if there are multiple good matches, test for a fork.
		# this is done with pairwise comparisons between matches - if do NOT align, there is a fork.
		nofork <- TRUE
		if (length(match_seqs) > 1){
			sub_match_inds <- 1:length(match_seqs)
			match_comparison_table <- unique_combos(x=sub_match_inds, y=sub_match_inds)
			good_align <- rep(F, nrow(match_comparison_table))
			for(j in 1:nrow(match_comparison_table)){
				seq1_ind <- match_comparison_table[[1]][j]
				seq2_ind <- match_comparison_table[[2]][j]
				good_align[j] <- testalign(match_seqs[seq1_ind], match_seqs[seq2_ind], search_word, return_s1_leftside=FALSE)
			}
			# if not all comparisons are good, there's a fork
			if( ! all(good_align)){nofork <- FALSE}
		}
		
		# if there's actually a match...
		if(sum(match_good) > 0){
			# add the first sequence, since it's the longest
			leftside_toadd <- testalign(match_seqs[1], seq_i, search_word, return_s1_leftside=TRUE)
			if(leftside_toadd == "FAIL"){stop("critical mission failure: somehow a sequence went from being GREAT, to being TERRIBLE. Buy new RAM?")}
			
			
			# update sequence that was added to
			if(direction == "right"){
				seqs_df$seq[i] <- rc(paste(leftside_toadd, seq_i, collapse=""))
			}else{
				seqs_df$seq[i] <- paste(leftside_toadd, seq_i, collapse="")
			}
		
			# update sequence data (name, length)
			seqs_df$length[i] <- nchar(seqs_df$seq[i])
			addednames_i <- paste(seqs_df$name[match_seqs_indices], collapse=",")
			seqs_df$name[i] <- paste(seqs_df$name[i], addednames_i, sep=",")
			
			#print("WIN!!")
			
			# If no fork, remove matched seq rows from seqs_df
			if(nofork){
				seqs_df <- seqs_df[-match_seqs_indices, ]
			# otherwise, just remove the one that was added
			}else{
				seqs_df <- seqs_df[-(match_seqs_indices[1]), ]
			}
					
		}
		
	
	}

	return(seqs_df)
}

simplify_contigs <- function(seqs_df){
	# for each contig, find and combine matches, then remove matches from db
	nseqs_before <- nrow(seqs_df)
	pb <- txtProgressBar(min = 0, max = nrow(seqs_df), style = 3)
	for(i in 1:(nseqs_before - 1)){
		
		if(i >= nrow(seqs_df)){break}
		
		seqs_df <- extend_contig(seqs_df, contig_index=i, word_length, direction="left")
		seqs_df <- extend_contig(seqs_df, contig_index=i, word_length, direction="right")
		
		# progress bar
		setTxtProgressBar(pb, i)
	}
	return(seqs_df)
}

write_seqs_df_to_fasta <- function(seqs_df, output_fp){
	for(i in 1:nrow(seqs_df)){
		write(paste(">", seqs_df$name[i]), file=output_fp, append=(i > 1) )
		write(seqs_df$seq[i], file=output_fp, append=TRUE)	
	}
}

