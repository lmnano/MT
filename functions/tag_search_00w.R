#tag and primer search in parallel
#making it faster by optimizing reverse complements
#further adapting for indels

tag_search_par = function(input = sequences_sample, tags = tags_uDNA, tag_mismatch = 2, search_length = 24, ncores = detectCores(), indels = TRUE){
  
  require(parallel)
  require(doParallel)
  require(foreach)
  require(Biostrings)
  
  #setting search range
  search_range = 1:search_length
  
  #searching for forward tags in parallel
  
  cl = makeCluster(ncores)
  clusterEvalQ(cl = cl, library(Biostrings))
  clusterExport(cl, c("tags", "input", "tag_mismatch", "search_range"), envir = environment())
  
  hits_tag_f = parSapply(cl = cl, input, function(X){
    y = X[search_range]
    countPDict(pdict = tags, subject = y, max.mismatch = tag_mismatch, with.indels = indels)
  })
  stopCluster(cl)
  
  registerDoParallel(cores = ncores)
  
  hits_tag_f_filter2 = foreach (i = 1:ncol(hits_tag_f), .combine = cbind, .packages = "Biostrings", .export = c("tags", "input", "tag_mismatch", "search_range")) %dopar% {
    #hits_tag_f_filter2 = foreach (i = 1:ncol(hits_tag_f), j = 1:nrow(hits_tag_f), .combine = cbind) %do% {
    
    filter_vector = rep(0, nrow(hits_tag_f))
    
    for (j in 1:nrow(hits_tag_f)) {
      # if (hits_tag_f[j,i] == 0) {
      #   filter_vector[[j]] = 0
      # }
      #this checks if found tag is in the correct place in the sequence
      if (hits_tag_f[[j,i]] == 1) {
        tag_filter = matchPattern(pattern = tags[[j]], subject = input[[i]][search_range], max.mismatch = tag_mismatch, with.indels = indels)
        tag_filter = as.data.frame(tag_filter@ranges)
        #if it is not it sets value in the matrix to 0
        if (indels == TRUE) {
          if (tag_filter$start %in% c(3,4,5) == FALSE) {
            filter_vector[[j]] = 0
          }else{
            filter_vector[[j]] = 1}
        }else{
          if (tag_filter$start != 4) {
            filter_vector[[j]] = 0
          }else{filter_vector[[j]] = 1}
        }
        
      }
      #in places where the same tag was found more than once in a sequence this checks if exactly one of the hits is in the correct position, if yes it sets number of found tags to 1 if not to 0
      #theoretical option where more than one hit is in the correct position also equals to 0
      if (hits_tag_f[[j,i]] > 1) {
        tag_filter = matchPattern(pattern = tags[[j]], subject = input[[i]][search_range], max.mismatch = tag_mismatch, with.indels = indels)
        tag_filter = as.data.frame(tag_filter@ranges)
        if (indels == TRUE) {
          if (sum(tag_filter$start %in% c(3,4,5)) == 1) {
            filter_vector[[j]] = 1
          }else{
            filter_vector[[j]] = 0}
        }else{
          if (sum(tag_filter$start == 4) == 1) {
            filter_vector[[j]] = 1
          }else{
            filter_vector[[j]] = 0}
        }
        
      }
    }
    filter_vector
  }
  
  registerDoSEQ()
  
  #searching for reverse tags in parallel
  
  tagsRC = reverseComplement(tags)
  
  cl = makeCluster(ncores)
  clusterEvalQ(cl = cl, library(Biostrings))
  clusterExport(cl, c("tags", "input", "tag_mismatch", "search_range", "tagsRC"), envir = environment())
  
  
  hits_tag_r = parSapply(cl = cl, input, function(X){
    y = X[(length(X)-search_length):length(X)]
    countPDict(pdict = tagsRC, subject = y, max.mismatch = tag_mismatch, with.indels = indels)
  })
  
  stopCluster(cl)
  
  registerDoParallel(cores = ncores)
  
  hits_tag_r_filter2 = foreach (i = 1:ncol(hits_tag_r), .combine = cbind, .packages = "Biostrings", .export = c("tags", "input", "tag_mismatch", "search_range", "tagsRC")) %dopar% {
    filter_vector = rep(0, nrow(hits_tag_f))
    
    for (j in 1:nrow(hits_tag_r)) {
      #this checks if found tag is in the correct place in the sequence...
      if (hits_tag_r[j,i] == 1) {
        tag_filter = matchPattern(pattern = tagsRC[[j]], subject = input[[i]][(width(input[i])-search_length):width(input[i])], max.mismatch = tag_mismatch, with.indels = indels)
        tag_filter = as.data.frame(tag_filter@ranges)
        #...if it is not it sets value in the matrix to 0
        if (indels == TRUE) {
          if (tag_filter$end %in% c(search_length - 2, search_length - 1, search_length - 3) == FALSE) {
            filter_vector[[j]] = 0
          }else{
            filter_vector[[j]] = 1}
        }else{
          if (tag_filter$end != search_length - 2) {
            filter_vector[[j]] = 0
          }else{
            filter_vector[[j]] = 1}
        }
        
      }
      #in places where the same tag was found more than once in a sequence this checks if exactly one of the hits is in the correct position, if yes it sets number of found tags to 1 if not to 0
      #theoretical option where more than one hit is in the correct position also equals to 0
      if (hits_tag_r[j,i] > 1) {
        tag_filter = matchPattern(pattern = tagsRC[[j]], subject = input[[i]][(width(input[i])-search_length):width(input[i])], max.mismatch = tag_mismatch, with.indels = indels)
        tag_filter = as.data.frame(tag_filter@ranges)
        if (indels == TRUE) {
          if (sum(tag_filter$end %in% c(search_length - 2, search_length - 1, search_length - 3)) == 1) {
            filter_vector[[j]] = 1
          }else{
            filter_vector[[j]] = 0}
        }else{
          if (sum(tag_filter$end == search_length - 2) == 1) {
            filter_vector[[j]] = 1
          }else{
            filter_vector[[j]] = 0}
        }
        
      }
    }
    filter_vector
  }
  registerDoSEQ()
  
  hits_tag_sum = hits_tag_f_filter2 + hits_tag_r_filter2
  
  hits_tag_f <<- hits_tag_f
  hits_tag_f_filter <<- hits_tag_f_filter2
  hits_tag_r <<- hits_tag_r
  hits_tag_r_filter <<- hits_tag_r_filter2
  hits_tag_sum <<- hits_tag_sum
  
}

primer_search_par = function(input = sequences_sample, primers = primers_uDNA, primer_mismatch = 3, search_length = 50, ncores = detectCores(), indels = TRUE){
  
  require(parallel)
  require(Biostrings)
  
  cl = makeCluster(ncores)
  clusterEvalQ(cl = cl, library(Biostrings))
  clusterExport(cl, c("tags", "input", "primer_mismatch", "search_length", "primers"), envir = environment())
  
  
  hits_primer_f2 = parSapply(cl = cl, input, function(X){
    y = X[1:search_length]
    countPDict(pdict = primers, subject = y, max.mismatch = primer_mismatch, with.indels = indels)
  })
  
  primersRC = reverseComplement(primers)
  
  hits_primer_r2 = parSapply(cl = cl, input, function(X){
    y = X[(length(X)-search_length):length(X)]
    countPDict(pdict = primersRC, subject = y, max.mismatch = primer_mismatch, with.indels = indels)
  })
  
  stopCluster(cl)
  
  
  hits_primer_sum = hits_primer_f2 + hits_primer_r2
  
  
  hits_primer_f <<- hits_primer_f2
  hits_primer_r <<- hits_primer_r2
  hits_primer_sum <<- hits_primer_sum
  
}  
#tag_search_par()