#Error simulation for DNA reads nicely wrapped
#This version uses normal strings for simulation, not DNAstrings and works only on one input .fastq file
#Does indels now
#Still needs some optimization
#Has setting for ratio between substitutions, insertions and deletions

error_sim = function(x = "./data/demo_data/", file = "DAB053demo01.tag16.fastq", y = "./data/error_data/", err_rate = 0.1, substitutions = 0.34, insertions = 0.33, deletions = 0.33) {
  
  stopifnot(substitutions + insertions + deletions == 1)
  stopifnot(err_rate < 1 & err_rate > 0)
  
  require(Biostrings)
  require(ShortRead)
  require(stringr)
  
  #data input
  reads = readFastq(dirPath = x, pattern = file)
  sequences = sread(reads)
  
  #Splitting the sequences into separate letters
  # vseq = as.vector(sequences)
  # vlseq= as.list(vseq)
  # vlseq_split = sapply(vlseq, str_split, pattern = "")
  
  #Same thing in one line
  vlseq_split = sapply(as.list(as.vector(sequences)), str_split, pattern = "")
  
  #bases to be used as replacements
  #base = c("A", "C", "T", "G", "-", "D")
  
  base = c(rep(c("A", "C", "T", "G"), times = (100*substitutions)), rep(rep("-", 4), times = (100*deletions)), rep(rep("D", 4), times = (100*insertions)))
  
  #Error simulation
  errdf = vlseq_split
  
  for (j in 1:length(vlseq_split)) {
    
    n_base_repl = round(err_rate * length(vlseq_split[[j]]))
    #n_base_repl = round(n_base_repl_full * (5/6))
    
    #if error rate is too low and/or dna reads are too short and as.integer function makes 0 errors, this makes 1 error
    if (n_base_repl == 0){
      n_base_repl = 1
    }
    
    #choosing random base number 
    rnd_n = sample(c(1:length(vlseq_split[[j]])), size = n_base_repl)
    #rnd_n = sample(rnd_n_full, size = round((5/6)*length(rnd_n_full)))
    
    rnd_b = sample(base, n_base_repl, replace = TRUE)
    
    #checking if switched bases are the same as original, if yes, ne next base in base vector is chosen, for the last base (G), the fitst one (A) is chosen
    
    #message(sprintf("Processing sequence %d/%d", j, length(vlseq_split)))
    
    for (i in 1:n_base_repl) {
      if ((vlseq_split[[j]][rnd_n[i]] == rnd_b[i]) & (rnd_b[i] != "G")){
        pos_base = match(rnd_b[i], base)
        rnd_b[i] = base[pos_base + 1]
      }else if((vlseq_split[[j]][rnd_n[i]] == rnd_b[i]) & (rnd_b[i] == "G")){
        rnd_b[i] = "A"
      }else if (rnd_b[i] == "D") {
        insert = sample(base[1:4], size = 1)
        rnd_b[i] = str_c(vlseq_split[[j]][rnd_n[i]], insert)
      }
    }
    
    # n_base_insert = setdiff(n_base_repl_full, n_base_repl)
    # rnd_n_insert = setdiff(rnd_n_full, rnd_n)
    # 
    # for (k in 1:n_base_insert) {
    #   rnd_insert[k] = str_c(vlseq_split[[j]][rnd_n_insert[k]], "a")
    # }
    
    errdf[[j]] = (replace(vlseq_split[[j]], rnd_n, rnd_b))
    
    
    
    errdf[[j]]
  }
  
  #making a DNAStringSet out of data with errors
  err_seq = lapply(errdf, base::paste, collapse='')
  err_seq = lapply(err_seq, str_remove_all, pattern = "-")
  errdf_len = lapply(err_seq, nchar)
  err_seq = lapply(err_seq, DNAString)
  err_seq = DNAStringSet(err_seq)
  
  #making BStringset for fastq quality reads
  #errdf_len = lapply(err_seq, length)
  errdf_k = as.list(rep("k", times = length(errdf)))
  
  for (i in 1:length(errdf)) {
    errdf_k[[i]] = rep("k", times = errdf_len[[i]])
    errdf_k
  }
  
  errdf_k = lapply(errdf_k, base::paste, collapse='')
  errdf_k = lapply(errdf_k, BString)
  errdf_k = BStringSet(errdf_k)
  errdf_k = SFastqQuality(errdf_k)
  
  
  #making a ShortReadq object for saving data as .fastq
  reads_err = ShortReadQ(sread = err_seq, quality = errdf_k, id = reads@id)
  
  #forming output file string
  out_file = str_replace(string = file, pattern = "(.*)\\.fastq", replacement = "\\1")
  out_err = paste("error", as.character(err_rate*100), sep = "")
  out_sub = paste("s", as.character(substitutions*100), sep = "")
  out_ins = paste("i", as.character(insertions*100), sep = "")
  out_del = paste("d", as.character(deletions*100), sep = "")
  out = paste(out_file, out_err, "indel", out_sub, out_ins, out_del, "fastq", sep = ".")
  out_final = paste(y, out, sep = "")
  
  #saving data
  
  writeFastq(object = reads_err, file = out_final, mode = "w", compress = TRUE)
  #reads_err <<- reads_err
  return()
  
}

#system.time(error_sim())
