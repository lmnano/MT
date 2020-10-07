#Function for 8 base tag identification
#Requres other functions for dara input (import_functions.R, tag_seqrch_v .R) and some packages to work
#Parallel version with apply loops, FORK clusters


primer_tag_id_parF = function(ftag_hits = hits_tag_f_filter, rtag_hits = hits_tag_r_filter, fprimer_hits = hits_primer_f, rprimer_hits = hits_primer_r, tags = tags_u, ngsfilter = ngsfilter_import, primers_u1 = primers_u, sequences_sample1 = sequences_sample, ncores = detectCores(), indels = TRUE){
  
  require(parallel)
  require(dplyr)
  
  #First tag identification
  
  #Matching found forward tags with tag sequences.
  
  coorf = which(ftag_hits == 1, arr.ind = TRUE)
  cl = makeCluster(ncores, type = "PSOCK")
  clusterExport(cl, c("coorf", "tags"), envir = environment())
  
  match_tag_f2 = t(parApply(cl = cl, coorf, 1, function(X){
    c(X[2], unname(X[1]), tags[[X[1]]])
  }))
  
  stopCluster(cl)
  
  colnames(match_tag_f2) = c("sequence_num", "f_tag_num", "f_tag")
  coorr = which(rtag_hits == 1, arr.ind = TRUE)
  
  cl = makeCluster(ncores, type = "PSOCK")
  clusterExport(cl, c("coorr", "tags"), envir = environment())
  
  match_tag_r2 = t(parApply(cl = cl, coorr, 1, function(X){
    c(X[2], unname(X[1]), tags[[X[1]]])
  }))
  
  stopCluster(cl)
  
  #Matching reverse tags with tag sequences.
  
  colnames(match_tag_r2) = c("sequence_num", "r_tag_num", "r_tag")
  
  #Combining both data frames.
  
  match_tag2 = merge(as.data.frame(match_tag_f2, stringsAsFactors = FALSE), as.data.frame(match_tag_r2, stringsAsFactors = FALSE), sort = FALSE)
  
  match_tag2$tag_combo_num = paste(match_tag2$f_tag_num, match_tag2$r_tag_num, sep = ":")
  match_tag2$tag_combo = paste(match_tag2$f_tag, match_tag2$r_tag, sep = ":")
  
  #Primer identification
  
  #Forward primer first, same as tag
  
  coorpf = which(fprimer_hits == 1, arr.ind = TRUE)
  cl = makeCluster(ncores, type = "PSOCK")
  clusterExport(cl, c("coorpf", "primers_u1"), envir = environment())
  
  match_primer_f2 = t(parApply(cl = cl, coorpf, 1, function(X){
    c(X[2], unname(X[1]), primers_u1[[X[1]]])
  }))
  
  stopCluster(cl)
  
  colnames(match_primer_f2) = c("sequence_num", "f_primer_num", "f_primer")
  coorpr = which(rprimer_hits == 1, arr.ind = TRUE)
  
  #Reverse primer
  cl = makeCluster(ncores, type = "PSOCK")
  clusterExport(cl, c("coorpr", "primers_u1"), envir = environment())
  
  match_primer_r2 = t(parApply(cl = cl, coorpr, 1, function(X){
    c(X[2], unname(X[1]), primers_u1[[X[1]]])
  }))
  
  stopCluster(cl)
  
  colnames(match_primer_r2) = c("sequence_num", "r_primer_num", "r_primer")
  
  #Combining data
  
  #Combining sequence numbers with tags and primers in one data frame.
  
  primer_tag_f2 = na.omit(left_join(x = as.data.frame(match_primer_f2, stringsAsFactors = FALSE), y = match_tag2))
  primer_tag_r2 = na.omit(left_join(x = as.data.frame(match_primer_r2, stringsAsFactors = FALSE), y = match_tag2))
  
  primer_tag2 = full_join(x = primer_tag_f2, y = primer_tag_r2)
  
  primer_tag2  = subset(primer_tag2, select = c("sequence_num", "f_primer", "r_primer", "tag_combo"))
  
  primer_tag_ngs_f2 = left_join(x = primer_tag2, y = ngsfilter, by = c("tag_combo" = "V3", "f_primer" = "V4"))
  primer_tag_ngs_r2 = left_join(x = primer_tag2, y = ngsfilter, by = c("tag_combo" = "V3", "r_primer" = "V5"))
  
  primer_tag_ngs_f2 = primer_tag_ngs_f2[,1:6]
  primer_tag_ngs_r2 = primer_tag_ngs_r2[,1:6]
  
  primer_tag_ngs2 = left_join(primer_tag_ngs_f2, primer_tag_ngs_r2, by = c("sequence_num", "f_primer", "r_primer", "tag_combo"))
  
  primer_tag_ngs2$locus = coalesce(primer_tag_ngs2$V1.x, primer_tag_ngs2$V1.y)
  primer_tag_ngs2$sample = coalesce(primer_tag_ngs2$V2.x, primer_tag_ngs2$V2.y)
  
  primer_tag_ngs2 = subset(primer_tag_ngs2, select = c("sequence_num", "f_primer", "r_primer", "tag_combo", "locus", "sample"))
  primer_tag_ngs2$sequence_num = as.integer(primer_tag_ngs2$sequence_num)
  
  #Adding missing primers
  
  primer_tag_ngs2 = left_join(x = primer_tag_ngs2, y = primers, by = c("locus" = "locusname"))
  primer_tag_ngs2 = primer_tag_ngs2[,c(1,4,5,6,7,8)]
  
  primer_tag_ngs2 = na.omit(primer_tag_ngs2)
  
  #Adding actual sequences to the data frame
  #this assumes both tags are of the same length
  #cl = makeCluster(ncores, type = "PSOCK")
  
  if (indels == TRUE) {
    
    #primer_tag_ngs2$sequence = parApply(cl = cl, primer_tag_ngs2, 1, function(X){
    primer_tag_ngs2$sequence = apply(primer_tag_ngs2, 1, function(X){
        
      c(as.character(sequences_sample1[[as.integer(X[names(X) == "sequence_num"])]]))
    })
    
  }else{
    #primer_tag_ngs2$sequence = parApply(cl = cl, primer_tag_ngs2, 1, function(X){
    primer_tag_ngs2$sequence = apply(primer_tag_ngs2, 1, function(X){
        
      c(as.character(sequences_sample1[[as.integer(X[names(X) == "sequence_num"])]][
        (1 + 3 + (((nchar(X[names(X) == "tag_combo"])) - 1) / 2)
         + nchar(X[names(X) == "primer1"]))
        :(width(sequences_sample1[as.integer(X[names(X) == "sequence_num"])])
          - (3 + (((nchar(X[names(X) == "tag_combo"])) - 1) / 2))
          - nchar(X[names(X) == "primer2"]))
      ]))
    })
  }
  
  #stopCluster(cl)
  
  primer_tag_ngs2$seq_length = nchar(primer_tag_ngs2$sequence)
  
  #primer_tag_ngs3 <<- primer_tag_ngs2
  primer_tag_ngs <<- primer_tag_ngs2
  
  
  
}
primer_tag_id_parF()