#final pipeline in R script

#Data import
source("./functions/import_functions_00.R")
demo_data_import(subset_seq = F, subset_size = 100, input_folder = "./data/error_data/", input_file = "DAB053demo01.tag16.error2.indel.s34.i33.d33.NEW.fastq")
tag_import(input_path = "./data/tags_primers/UrsusNGSPrimersCrosbreeding.Ljubljana.18Jul2017_16tag.switch.xlsx", input_sheet = "Tag crossbreeding")
primer_import()
ngsfilter_import(input_path = "./data/ngsfilter/DAB053demo01.tag16.ngsfilter")


#Tag and primer search
if (.Platform$OS.type == "unix") {
  source("./functions/tag_search_00.R")
  
}else if (.Platform$OS.type == "windows"){
  source("./functions/tag_search_00w.R")
  
}

system.time(tag_search_par(indels = T))
system.time(primer_search_par(indels = T))

library(data.table)
output = sprintf("./data/tag_primer_hits/hits_primer_f_%s.txt", data_name)
fwrite(x = hits_primer_f, file = output)
output = sprintf("./data/tag_primer_hits/hits_primer_r_%s.txt", data_name)
fwrite(x = hits_primer_r, file = output)
output = sprintf("./data/tag_primer_hits/hits_primer_sum_%s.txt", data_name)
fwrite(x = hits_primer_sum, file = output)
output = sprintf("./data/tag_primer_hits/hits_tag_f_%s.txt", data_name)
fwrite(x = hits_tag_f, file = output)
output = sprintf("./data/tag_primer_hits/hits_tag_r_%s.txt", data_name)
fwrite(x = hits_tag_r, file = output)
output = sprintf("./data/tag_primer_hits/hits_tag_f_filter_%s.txt", data_name)
fwrite(x = hits_tag_f_filter, file = output)
output = sprintf("./data/tag_primer_hits/hits_tag_r_filter_%s.txt", data_name)
fwrite(x = hits_tag_r_filter, file = output)
output = sprintf("./data/tag_primer_hits/hits_tag_sum_%s.txt", data_name)
fwrite(x = hits_tag_sum, file = output)


#Primer and tag id
if (.Platform$OS.type == "unix") {
  source("./functions/primer_tag_id_00.R")
  
}else if (.Platform$OS.type == "windows"){
  source("./functions/primer_tag_id_00w.R")
  
}
system.time(primer_tag_id_parF(indels = T))

# library(data.table)
output = sprintf("./data/primer_tag_ngs/primer_tag_ngs_%s.txt", data_name)
fwrite(x = primer_tag_ngs, file = output)

#Alignment of sample-cluster groups and length correction
source("./functions/alignment_00.R")
system.time(alignment())

#source("./functions/alignment_00w.R")
#system.time(alignment())

#library(data.table)
output = sprintf("./data/fixSeq_alignment_output/fixSeq_%s.txt", data_name)
fwrite(x = fixSeq, file = output)


# Clustering and consesnsus
source("./functions/clustering_consensus_00.R")
system.time(clustering_par(indel = T, sorted_data = fixSeq))

output = sprintf("./data/grouped_data_consensus/grouped_data_consensus_%s.txt", data_name)
fwrite(x = grouped_data_consensus_par_mclapply, file = output)


#Fishbone allele calling
source("./functions/fishbone_prepare_00.R")
fishbone_prepare(filter = F)

library(fishbone)

out2 <- mclapply(fishbone_input, FUN = callAllele, tbase = mt, verbose = FALSE)#, clean = FALSE)


#fishbone_out <- do.call(rbind, out)
fishbone_out <- do.call(rbind, out2)

library(data.table)
output = sprintf("./data/fishbone_out/fishbone_out_%s.txt", data_name)
fwrite(x = fishbone_out, file = output)


