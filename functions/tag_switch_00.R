#script for switching tags in demo data
#tags need to be generated with appropriate tag generator for the script to import and isolate new tags properly
#default demo data for switching tags uses 8 bp tags
#default new tags are 16 bp tags used in later pipeline


tag_switch = function(new_tags = "./data/new_tags/barcode_16.txt", output = "./data/demo_data/", import_function_source = "./functions/import_functions_00.R", demo_data_input_folder = "./data/demo_data", demo_data_input_file = "DAB053_demo01.fastq", old_tag_input_path = "./data/tags_primers/Ursus\ NGS\ Primers\ list\ with\ tags.xlsx", old_tag_input_sheet = "Tag crossbreeding "){
  
  #importing required libraries and sourcing functions
  require(dplyr)
  require(Biostrings)
  require(ShortRead)
  require(readxl)
  require(stringr)
  
  source(file = import_function_source)
  
  
  #importing new tags
  tag_import = read.table(new_tags, skip = 1, stringsAsFactors = FALSE)
  
  #importing data
  demo_data_import(input_folder = demo_data_input_folder, input_file = demo_data_input_file, subset_seq = FALSE)
  
  #importing old tags
  tag_import(input_path = old_tag_input_path, input_sheet = old_tag_input_sheet)
  
  #sorting new tags
  tags_new = tag_import %>%
    filter(grepl(">adA2_.*", V1) == TRUE) %>%
    mutate(tags = gsub(">adA2_(.*)", "\\1", V1)) %>%
    select(tags) %>%
    mutate(tags = toupper(tags))
  
  tags_newDNA = sapply(tags_new$tags, DNAString)
  tags_newDNA = DNAStringSet(tags_newDNA)
  
  #making reverse complements and preparing old and new tags for switching
  seq_char = as.data.frame(as.character(sequences))
  
  tags_u12DNA = tags_newDNA[1:12]
  tags_u12DNArc = reverseComplement(tags_u12DNA)
  
  tags_uDNArc = reverseComplement(tags_uDNA)
  
  tags12PP1 = tags_new[1:12,]
  tags12PP1rc = unname(as.character(tags_u12DNArc))
  
  tags8PP1 = toupper(tags_u[,])
  tags8PP1rc = as.character(tags_uDNArc)
  
  #switching tags in sequences
  
  seq_ch = as.character(sequences)
  
  seq12 = seq_ch
  
  for (i in 1:length(seq12)) {
    for (j in 1:length(tags8PP1)) {
      seq12[i] = gsub(pattern = tags8PP1[j], replacement = tags12PP1[j], x = seq12[i])
    }
    if (i %% 10000 == 0) {
      message(sprintf("Sequence %d/%d.", i, length(sequences)))
    }
    
  }
  
  for (i in 1:length(seq12)) {
    for (j in 1:length(tags8PP1rc)) {
      seq12[i] = gsub(pattern = tags8PP1rc[j], replacement = tags12PP1rc[j], x = seq12[i])
    }
    if (i %% 10000 == 0) {
      message(sprintf("Reverse complement sequence %d/%d.", i, length(sequences)))
    }
  }
  
  seq12DNA = DNAStringSet(seq12)
  
  #making quality reads for writing in fastq
  
  len12 = unname(sapply(seq12, nchar))
  k = as.list(rep("k", times = length(seq12)))
  
  for (i in 1:length(seq12)) {
    k[[i]] = rep("k", times = len12[[i]])
    k
  }
  
  k = lapply(k, base::paste, collapse='')
  k = lapply(k, BString)
  k = BStringSet(k)
  k = SFastqQuality(k)
  
  #writing output file and path string
  out_tag_length = paste("tag", as.character(nchar(tags_new[1,1])), sep = "")
  out_file = str_replace(string = demo_data_input_file, pattern = "(.*)\\.fastq", replacement = "\\1")
  out = paste(out_file, out_tag_length, "fastq", sep = ".")
  out_final = paste(output, out, sep = "")
  
  #writing fastq
  reads_seq12 = ShortReadQ(sread = seq12DNA, quality = k, id = reads@id)
  writeFastq(object = reads_seq12, file = out_final, mode = "w", compress = FALSE)
  
}

#tag_switch()