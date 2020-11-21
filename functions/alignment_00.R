#alignmnet and length fixing function
#line 22 is for testing only

alignment = function(input = primer_tag_ngs){
  
  require(dplyr)
  require(Biostrings)
  require(stringr)
  require(parallel)
  
  #filtering out multiple matches and splittin in sample-locus groups
  # primer_tag_ngs2 = input %>%
  #   group_by(sequence) %>%
  #   filter(n() == 1) %>%
  #   ungroup() %>%
  #   group_by(sample, locus) %>%
  #   arrange(.by_group = TRUE) %>%
  #   add_count() %>%
  #   group_split()
 
  primer_tag_ngs2 = primer_tag_ngs %>%
    group_by(sequence) %>%
    group_by(sequence_num) %>%
    distinct() %>%
    add_count() %>%
    filter(n == 1) %>%
    ungroup() %>%
    select(-n) %>%
    group_by(sample, locus) %>%
    arrange(.by_group = TRUE) %>%
    add_count() %>%
    group_split()
  
  
  #for testing only
  #primer_tag_ngs2 = primer_tag_ngs2[1:3]
  
  #making locus03 fake reference sequence
  ref03flank1 = DNAString("AAATCCTGTAACAAAT")
  ref03flank2 = DNAString("CTC")
  
  motif03 = "CTAT"
  repeat03 = rep(motif03, times = 25)
  repeat03 = DNAString(paste(repeat03, collapse = ""))
  
  fictional03 = DNAString(c(ref03flank1, repeat03, ref03flank2))
  
  #making locus06 fake sequence
  #ref06flank1 = DNAString("")
  ref06flank2 = DNAString("AAAAGAAGACAGATTGTAA")
  
  motif06 = "AAGG"
  repeat06 = rep(motif06, times = 25)
  repeat06 = DNAString(paste( repeat06, collapse = ""))
  
  fictional06 = DNAString(c(repeat06, ref06flank2))
  
  #making locus17 flanking sequence
  
  ref17flank1 = DNAString("TAGCTGGTTTTCTTTTT")
  ref17flank2 = DNAString("GATGGATATTTATTTCT")
  
  motif17 = "CTTT"
  repeat17 = rep(motif17, times = 25)
  repeat17 = DNAString(paste(repeat17, collapse = ""))  
  
  fictional17 = DNAString(c(ref17flank1, repeat17, ref17flank2))
  
  start = DNAString("GGG")
  end = DNAString("TTT")
  
  fixSeq = mclapply(primer_tag_ngs2, function(y){
    
    #extracting tags and primers for appropriate reference sequence
    fprimer = DNAString(unique(y$primer1))
    rprimer = reverseComplement(DNAString(unique(y$primer2)))
    tag = unique(y$tag_combo)
    tag = str_split(tag, pattern = ":")
    tag1 = DNAString(tag[[1]][1])
    tag2 = reverseComplement(DNAString(tag[[1]][2]))
    
    if (y$locus[1] == "UA_MxRout1_03") {
      fictional = DNAStringSet(c(start, tag1, fprimer, fictional03, rprimer, tag2, end))
    }else if (y$locus[1] == "UA_MxRout1_06"){
      fictional = DNAStringSet(c(start, tag1, fprimer, fictional06, rprimer, tag2, end))
    }else if (y$locus[1] == "UA_MxRout1_17"){
      fictional = DNAStringSet(c(start, tag1, fprimer, fictional17, rprimer, tag2, end))
    }
    
    #names(fictional) = "reference"
    
    #extracting sequences
    err_seqs = DNAStringSet(y$sequence)
    
    fakeHead1 = fictional[[1]][1:round(width(fictional)/2)]
    fakeTail1 = fictional[[1]][round(width(fictional)/2):width(fictional)]
    
    #alignment and sequence fixing
    errseqFix_all = sapply(err_seqs, function(x, fakeHead, fakeTail){
      
      
      
      #pairwise alignments
      alHead = pairwiseAlignment(x, fakeHead, type = "local", gapOpening = 5)
      
      alTail = pairwiseAlignment(x, fakeTail, type = "local", gapOpening = 5)
      
      #marking bases to be deleted from the error sequence
      gap_matchesH = matchPattern(pattern = "-", subject = (DNAStringSet(alHead@subject))[[1]])
      gap_matchesH = gap_matchesH@ranges@start
      repl_letters = rep("x", times = length(gap_matchesH))
      delHead =  length(gap_matchesH)
      
      errseqHead = as.character(alHead@pattern)
      errseqHead = str_split(errseqHead, pattern = "")
      errseqHead = replace(errseqHead[[1]], gap_matchesH, repl_letters)
      errseqHead = paste(errseqHead, collapse = "")
      #errseqHead = str_remove_all(errseqHead, pattern = "x")
      
      gap_matchesT = matchPattern(pattern = "-", subject = (DNAStringSet(alTail@subject))[[1]])
      gap_matchesT = gap_matchesT@ranges@start
      repl_letters = rep("x", times = length(gap_matchesT))
      delTail =  length(gap_matchesT)
      
      errseqTail = as.character(alTail@pattern)
      errseqTail = str_split(errseqTail, pattern = "")
      errseqTail = replace(errseqTail[[1]], gap_matchesT, repl_letters)
      errseqTail = paste(errseqTail, collapse = "")
      #errseqTail = str_remove_all(errseqTail, pattern = "x")
      
      #merging overlaping error sequence alignments
      errHeadend = alHead@pattern@range@width
      errTailstart = alTail@pattern@range@start
      errTailend = alTail@pattern@range@width
      
      head_clipped = substr(errseqHead, 1, errHeadend - abs(errTailstart - errHeadend))
      
      errseqFix = paste(head_clipped, errseqTail, sep ="")
      
      #removing bases marked to be deleted
      errseqFix = str_remove_all(errseqFix, pattern = "x")
      
      errseqFix
      
      
    }, fakeHead = fakeHead1, fakeTail = fakeTail1)
    
    
    fix_seq_length = nchar(errseqFix_all)
    lenFixseq = as.data.frame(errseqFix_all, stringsAsFactors = FALSE)
    
    lenFixseq = lenFixseq %>%
      mutate(length = fix_seq_length) %>%
      mutate(len_fix = 1) %>%
      group_by(length) %>%
      add_count() %>%
      arrange(desc(n))
    
    base_length = lenFixseq$length[1]
    
    lenFixseq = group_split(lenFixseq)
    
    lenFixseq = lapply(lenFixseq, function(x){
      xlen = x$length[1]
      
      if (abs(xlen-base_length) %% 4 == 0) {
        x$len_fix = xlen
      }else if ((abs(xlen + 1 - base_length)) %% 4 == 0){
        x$len_fix = xlen + 1 
        #}else if ((abs(xlen + 2 - base_length)) %% 4 == 0){
        #  x$len_fix = xlen + 2
      }else if ((abs(xlen - 1 - base_length)) %% 4 == 0){
        x$len_fix = xlen - 1 
      }
      x
      
    })
    
    lenFixseq = do.call(rbind, lenFixseq)
    
    lenFixseq = lenFixseq %>%
      select(-n) %>%
      ungroup() %>%
      group_by(len_fix) %>%
      add_count()
    
    primer_tag_ngs2fix = y
    
    primer_tag_ngs2fix$sequence = lenFixseq$errseqFix_all
    primer_tag_ngs2fix$seq_length = lenFixseq$len_fix
    primer_tag_ngs2fix$n = lenFixseq$n
    primer_tag_ngs2fix = filter(primer_tag_ngs2fix, (n/nrow(primer_tag_ngs2fix)) > 0.1)
    primer_tag_ngs2fix = filter(primer_tag_ngs2fix, seq_length != 1)
    
    primer_tag_ngs2fix
    
    
  })
  
  
  #length value fixing
  #inserting length corrections for later grouping
  fixSeq = do.call(rbind, fixSeq)
  fixSeq <<- fixSeq
}

#system.time(alignment())