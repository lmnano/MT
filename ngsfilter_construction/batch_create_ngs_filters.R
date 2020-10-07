# This script will take aliquot plates, forward and reverse primer sequences, sequences of tags and create
# a combination
# Glossary:
# PP = primer plate, combination of tags and primers. Tags are unique for each well. Eight plates per library.
# AP = aliquot plate, each well holds its own sample.

# Make sure you change locations of files to read according to the plates you're using.

library(readxl)

project <- "ngsfilter_construction"
dir.AP <- sprintf("./%s/0_prep_ngsfilters/aliquot_plates", project) # folder where aliquote plates are located
dir.output <- sprintf("./%s/1_ngsfilters", project) # no trailing slash, where files are to be stored

# This file contains data which maps which aliquot plate comes from which library.
## I've changed platenemaes for demo data
platenames <- "190916_PlateNames_GATC_demo.xlsx"
map.AP <- sprintf("./%s/0_prep_ngsfilters/%s", project, platenames)
map.sheet <- "PCR_Plates_4Reps_2Reps"

# Combination of tags to determine sample position/identity.
# PP has columns position, slo, PP1, PP2, ... PP8 which designates which position (1-96) holds which forward
# and reverse tag combination.
combo.PP <- sprintf("./%s/0_prep_ngsfilters/UrsusNGSPrimersCrosbreeding.Ljubljana.18Jul2017_16tag.switch.xlsx", project)

if (!dir.exists(dir.output)) {
  stop(sprintf("The output folder %s you specified does not exists. 
               Make sure it exists beforehand running this script.", dir.output))
}

if (!dir.exists(dir.AP)) {
  stop(sprintf("The folder which should hold aliquot plates (%s) does't 
       appear to be there. Please check your paths.", dir.AP))
}

# 1. Load PP and AP data
# Primer names and forward/reverse sequences.
primers <- as.data.frame(read_excel(sprintf("./%s/0_prep_ngsfilters/input_primers_tags.xlsx", project), sheet = "primers"))

PP <- as.data.frame(read_excel(combo.PP, 
                               sheet = "Tag crossbreeding", skip = 43))

if (nrow(PP) == 0) {
  stop("It doesn't appear you've imported any primer plates. Please make sure the files are in place.")
}

# AP will hold links to files to aliquot plates. Data is arranged in columns. Each row holds its own
# sample (name) which will be used to construct .ngsfilter.
AP <- data.frame(location = list.files(dir.AP, pattern = "_[AB]\\d+\\.xls$", full.names = TRUE))
AP$name <- gsub("^.*([AB]\\d+)\\.xls", "\\1", AP$location)

if (nrow(AP) == 0) {
  stop("You have imported primer plates, but filtering failed. Make sure the regex expression in the above lines is correct.")
}

# 2. Find which AP is added to which PP.
pa.loc <- as.data.frame(read_excel(map.AP, sheet = map.sheet))

# select which libraries you wish to run through this script
# pa.loc <- droplevels(pa.loc[pa.loc$Library_BC %in% sprintf("DAB%02d", 13:24, sep = ""), ])

pa.loc <- split(pa.loc, f = pa.loc$Library_BC)

# For library, split by plate...
out <- sapply(pa.loc, FUN = function(x, PP, AP, primers, outloc = dir.output) {
  x.split <- split(x, f = 1:nrow(x))
  
  # ... and for each plate, construct NGS filter
  out <- sapply(x.split, FUN = function(y, PP, AP, primers) {
    # First find primers for PP in library
    find.pp <- y[, "Primer Plate"]
    PP.x <- na.omit(PP[, names(PP) %in% find.pp, drop = FALSE]) # columns plate and tagcombo
    stopifnot(ncol(PP.x) == 1)
    
    find.ap <- y[, "Aliquot Plate"]
    AP.x <- AP[AP$name %in% find.ap, ] # columns location and name
    AP.x <- read_excel(as.character(AP.x$location))
    
    pos.number <- sprintf("%03d", rep(1:nrow(AP.x), each = nrow(primers)))
    sample.pos.plate <- rep(AP.x$SPositionBC, each = nrow(primers))
    sample.pos.plate <- paste(sample.pos.plate, pos.number, find.pp, sep = "_")
    
    locus <- rep(primers$locusname, times = nrow(AP.x))
    
    tagcombo <- rep(PP.x[, 1], each = nrow(primers))
    
    primer1 <- rep(primers$primer1, times = nrow(PP.x))
    primer2 <- rep(primers$primer2, times = nrow(PP.x))
    
    out <- data.frame(locus, sample.pos.plate, tagcombo, primer1, primer2, ef = "F", at = "@")
  }, PP = PP, AP = AP, primers = primers, simplify = FALSE)
  
  out <- do.call(rbind, out)
  
  write.table(out, file = sprintf("%s/%s.ngsfilter", outloc, unique(x$Library_BC)), row.names = FALSE,
              col.names = FALSE, quote = FALSE, sep = "\t", fileEncoding = "UTF-8")
  out
}, PP = PP, AP = AP, primers = primers, outloc = dir.output, simplify = FALSE)

