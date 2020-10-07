#This is a clustering function for grouped sequences
#works with apply, alignment is uses all cores by default
#sorts clusters by size
#it works!!!

clustering_par = function(sorted_data = fixSeq, dist_matrix_method = "lv", cluster_n = 2, min_seq_num = 20, clustering = TRUE, indel = TRUE){
  
  if (cluster_n != 2) {
    stop("Only works on two clusters for now")
  }
  
  require(dplyr)
  require(Biostrings)
  require(DECIPHER)
  require(stringdist)
  
  #grouping the data
  grouped_data = sorted_data %>%
    group_by(sample, locus, seq_length) %>%
    arrange(.by_group = TRUE) %>%
    add_count() %>%
    filter(n >= cluster_n) %>%
    group_split()
  
  grouped_data_consensus = lapply(grouped_data, function(x){
    
    if (clustering == TRUE) {
      hicount_seq = x$sequence
      hicount_seq_len = length(hicount_seq)
      
      if (indel == TRUE) {
        hicount_seqDNA = RemoveGaps(DNAStringSet(hicount_seq))
        hicount_seqDNA = AlignSeqs(hicount_seqDNA, verbose = FALSE, processors = NULL)
        
        hicount_seq = as.character(hicount_seqDNA)
      }
      
      
      
      if (hicount_seq_len < min_seq_num) {
        hicount_cl1cut = rep(1, times = hicount_seq_len)
      }else{
        #making a distance matrix
        #message(sprintf("Making distance matrix for row %d of %d.", i, nrow(grouped_data)))
        hicount_sd1 = stringdistmatrix(hicount_seq, method = dist_matrix_method)
        
        #clustering step
        hicount_cl1 = hclust(hicount_sd1, method = "complete")
        hicount_cl1cut = cutree(hicount_cl1, k = cluster_n)
        
        clSum = sapply(1:cluster_n, function(x){
          out = sum(hicount_cl1cut == x)
          out
        })
        #if the clusters are too small they get merged
        if (((clSum[1]/hicount_seq_len) < 0.1) || ((clSum[2]/hicount_seq_len) < 0.1)) {
          hicount_cl1cut = rep(1, times = hicount_seq_len)
        }
      }
      
      locusf = x$locus[[1]]
      samplef = x$sample[[1]]
      seq_lenf = x$seq_length[[1]]
      nf = x$n[[1]]
      tagcombof = x$tag_combo[[1]]
      
      x$cluster = hicount_cl1cut
      
      cons_seqs = data.frame(sample = samplef, tag_combo = tagcombof, locus = locusf, seq_length = seq_lenf, n = nf, cluster = c(1:cluster_n), num_of_seqs_in_cluster = NA, consensus_seq = NA)
      
      for (k in 1:cluster_n) {
        cons_seqs[k, "num_of_seqs_in_cluster"] = sum(hicount_cl1cut == k)
        
        hicount_group_cl1 = filter(x, cluster == k)
        hicount_group_cl1 = DNAStringSet(hicount_group_cl1$sequence)
        
        hicount_cons1 = ConsensusSequence(hicount_group_cl1, threshold = 0.49)
        cons_seqs[k, "consensus_seq"] = as.character(hicount_cons1)
      }
      cons_seqs = filter(cons_seqs, num_of_seqs_in_cluster != 0)
      
      cons_seqs
    }else{
      locusf = x$locus[[1]]
      samplef = x$sample[[1]]
      seq_lenf = x$seq_length[[1]]
      nf = x$n[[1]]
      tagcombof = x$tag_combo[[1]]
      

      cons_seqs = data.frame(sample = samplef, tag_combo = tagcombof, locus = locusf, seq_length = seq_lenf, n = nf, consensus_seq = NA)
      

      hicount_group_cl1 = DNAStringSet(x$sequence)
      
      hicount_cons1 = ConsensusSequence(hicount_group_cl1, threshold = 0.49)
      cons_seqs[, "consensus_seq"] = as.character(hicount_cons1)

      cons_seqs
    }
    

  })
  grouped_data_consensus = do.call("rbind", grouped_data_consensus)
  grouped_data_consensus_par_mclapply <<- grouped_data_consensus
}

#clustering01_par_mclapply()