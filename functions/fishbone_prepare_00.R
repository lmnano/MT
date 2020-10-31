#Fishbone prepare function
#This function prepares data for allele calling with fishbone

fishbone_prepare = function(input = grouped_data_consensus_par_mclapply, library_path = "./ngsfilter_construction/0_prep_ngsfilters/190916_PlateNames_GATC_demo.xlsx", clustered_data = TRUE, filter = FALSE){
  
  require(dplyr)
  require(readxl)
  
  NGSlibrary = read_excel(path = library_path, sheet = "PCR_Plates_4Reps_2Reps")
  NGSlibrary = NGSlibrary$Library_BC
  
  if (clustered_data == TRUE) {
    if (filter == FALSE) {
      grouped_data_consensus_fb = input %>%
        #filter(cluster == 1) %>%
        mutate(locus = as.character(locus)) %>%
        mutate(locus = gsub(pattern = ".*(\\d{2})$", replacement = "\\1", x = locus)) %>%
        dplyr::rename(Marker = locus) %>%
        mutate(Plate = gsub(pattern = ".*\\w{2}(\\d{+})$", replacement = "\\1", x = sample)) %>%
        mutate(Position = gsub(pattern = ".*_(\\d*)_.*", replacement = "\\1", x = sample)) %>%
        mutate(Sample_Name = gsub(pattern = "^(.*)_.*_.*$", replacement = "\\1", x = sample)) %>%
        dplyr::rename(Read_Count = num_of_seqs_in_cluster) %>%
        #next line for NOcluster only, otherwise use line above
        #dplyr::rename(Read_Count = n) %>%
        dplyr::rename(length = seq_length) %>%
        dplyr::rename(Sequence = consensus_seq) %>%
        mutate(Sequence = tolower(as.character(Sequence))) %>%
        select(-c(sample, n, cluster)) %>%
        #next line for NOcluster only, otherwise use line above
        #select(-c(sample)) %>%
        
        #manually added needs to be automated
        mutate(Run_Name = NGSlibrary) %>%
        dplyr::rename(TagCombo = tag_combo)
      
    }else{
      grouped_data_consensus_fb = input %>%
        filter(cluster == 1) %>%
        mutate(locus = as.character(locus)) %>%
        mutate(locus = gsub(pattern = ".*(\\d{2})$", replacement = "\\1", x = locus)) %>%
        dplyr::rename(Marker = locus) %>%
        mutate(Plate = gsub(pattern = ".*\\w{2}(\\d{+})$", replacement = "\\1", x = sample)) %>%
        mutate(Position = gsub(pattern = ".*_(\\d*)_.*", replacement = "\\1", x = sample)) %>%
        mutate(Sample_Name = gsub(pattern = "^(.*)_.*_.*$", replacement = "\\1", x = sample)) %>%
        dplyr::rename(Read_Count = num_of_seqs_in_cluster) %>%
        #next line for NOcluster only, otherwise use line above
        #dplyr::rename(Read_Count = n) %>%
        dplyr::rename(length = seq_length) %>%
        dplyr::rename(Sequence = consensus_seq) %>%
        mutate(Sequence = tolower(as.character(Sequence))) %>%
        select(-c(sample, n, cluster)) %>%
        #next line for NOcluster only, otherwise use line above
        #select(-c(sample)) %>%
        
        #manually added needs to be automated
        mutate(Run_Name = NGSlibrary) %>%
        dplyr::rename(TagCombo = tag_combo)
      
    }
    
    
  }else{
    grouped_data_consensus_fb = input %>%
      mutate(locus = as.character(locus)) %>%
      mutate(locus = gsub(pattern = ".*(\\d{2})$", replacement = "\\1", x = locus)) %>%
      dplyr::rename(Marker = locus) %>%
      mutate(Plate = gsub(pattern = ".*\\w{2}(\\d{+})$", replacement = "\\1", x = sample)) %>%
      mutate(Position = gsub(pattern = ".*_(\\d*)_.*", replacement = "\\1", x = sample)) %>%
      mutate(Sample_Name = gsub(pattern = "^(.*)_.*_.*$", replacement = "\\1", x = sample)) %>%
      dplyr::rename(Read_Count = n) %>%
      dplyr::rename(length = seq_length) %>%
      dplyr::rename(Sequence = consensus_seq) %>%
      mutate(Sequence = tolower(as.character(Sequence))) %>%
      select(-c(sample)) %>%
      
      #manually added needs to be automated
      mutate(Run_Name = NGSlibrary) %>%
      dplyr::rename(TagCombo = tag_combo)
    
  }
  
  xz <- split(grouped_data_consensus_fb, f = list(grouped_data_consensus_fb$Sample_Name, grouped_data_consensus_fb$Marker, grouped_data_consensus_fb$Plate, grouped_data_consensus_fb$Run_Name, grouped_data_consensus_fb$Position))
  xz <-  xz[lapply(xz, nrow)>0]
  
  grouped_data_consensus_fb <<- grouped_data_consensus_fb
  fishbone_input <<- xz
}
#fishbone_prepare()
