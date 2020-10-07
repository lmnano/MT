#Functions for my demo data, primer and tag import
#They are written specifically for my demo data and other data I often use
#These functions requre certain libraries inportet to work


#import 8 tag demo data with errors
demo_data_import = function(input_folder = "./data/demo_data", input_file = "DAB053demo01.tag16.fastq", subset_seq = FALSE, subset_size = 100, filename = TRUE){
  require(ShortRead)
  
  reads = readFastq(dirPath = input_folder, pattern = input_file)
  sequences = sread(reads)
  
  if (subset_seq == TRUE) {
    sequences_sample = sample(sequences, subset_size)
  }else{
    sequences_sample = sequences
  }
  
  if (filename == TRUE) {
    data_name = sub(pattern = "(.*)\\.fastq", replacement = "\\1", x = input_file)
    data_name  <<- data_name
  }
  
  reads <<- reads
  sequences <<- sequences
  sequences_sample <<- sequences_sample
}

#import 8 base tags
tag_import = function(input_path = "./data/tags_primers/UrsusNGSPrimersCrosbreeding.Ljubljana.18Jul2017_16tag.switch.xlsx", input_sheet = "Tag crossbreeding", skip_rows = 43){
  require(readxl)
  require(Biostrings)
  require(stringr)
  
  #import tags
  tags = as.data.frame(read_excel(path = input_path, sheet = input_sheet, skip = skip_rows))
  
  #selecting tags for demo data
  tags_u = tags[,3]
  tags_u = matrix(str_split_fixed(as.matrix(tags_u), ":", 2), ncol = 1)
  tags_u = unique(tags_u)
  
  tags_uDNA = DNAStringSet(tags_u)
  
  tags <<- tags
  tags_u <<- tags_u
  tags_uDNA <<- tags_uDNA
}

#import primers

primer_import = function(input_path = "./data/tags_primers/input_primers_tags.xlsx", primer_sheet = "primers"){
  require(readxl)
  require(Biostrings)
  
  #import primers
  primers = as.data.frame(read_excel(input_path, sheet = primer_sheet))
  
  primers_u = matrix(unlist(primers[,2:3], use.names = FALSE))
  primers_uDNA = DNAStringSet(primers_u)
  
  primers <<- primers
  primers_u <<- primers_u
  primers_uDNA <<- primers_uDNA
}

ngsfilter_import = function(input_path = "./data/ngsfilter/DAB053demo01.tag16.ngsfilter"){
  ngsfilter8 = read.table(input_path, header = FALSE, stringsAsFactors = FALSE)
  
  ngsfilter_import <<- ngsfilter8
}


