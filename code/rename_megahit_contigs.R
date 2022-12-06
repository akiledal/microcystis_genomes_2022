#!/usr/bin/env Rscript

'Rename MegaHit contigs & export information from header as a .tsv file

Usage:
  rename_megahit_contigs.R <contigs> [options]

Options:
  -i CONTIGS --contigs=CONTIGS Uncompressed contigs .fa produced by MegaHit
  -h --help                   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
#arguments <- docopt(doc, args = c(" ~/GLAMR/data/omics/metagenomes/0308a1e0034cfa73ce078604e79e9ee7/assembly/megahit/final.contigs.fa"))
#print(arguments)

###########################
### Script starts here ####
###########################

#contigs <- "data/omics/metagenomes/4f5d4a3bdeb0797088688cfd48f1829c/assembly/megahit/final.contigs.fa"

library(tidyverse,quietly = TRUE)

if (!is.null(arguments$`--contigs`)) {
  contig_path <- arguments$`--contigs`
} else if(!is.null(arguments$contigs)){
  contig_path <- arguments$contigs
}

rename_contigs <- function(contigs){
  
  sample <- str_remove(contigs, ".*prodigal_mag/") %>% str_remove("\\.faa") %>% str_remove_all("_")
  contig_info_fp <- str_remove(contigs, "\\.faa") %>% paste0(.,"_contigs_info.tsv")
  
  contigs_df <- Biostrings::readAAStringSet(contigs) %>% 
    data.frame(header = names(.),
               seq = .)
  getwd()
  
  if(contigs_df$header[1] == paste0(sample,"_","1")){
    print("Contigs have already been processed")
    quit(save = "no", status = 1)
  }
  
  contigs_df <- contigs_df %>% 
    separate(header, into = c("orig_id","flag","approx_cov", "length"),sep = " ") %>% 
    mutate(across(c(flag,approx_cov, length), ~str_remove(.,".*="))) %>% 
    type_convert() %>%
    mutate(sample = sample,
           contig_id = paste0(sample,"_",row_number())) %>% 
    relocate(contig_id)
  
  contig_info <- contigs_df %>% 
    select(contig_id, orig_id, flag, approx_cov, length) %>% 
    write_tsv(contig_info_fp)

  out_path <- paste0(str_remove(contigs,"\\.faa"),".renamed.faa")
  
  for_export <- Biostrings::AAStringSet(contigs_df$seq)
  names(for_export) <- contigs_df$contig_id
  
  Biostrings::writeXStringSet(for_export,out_path)
}

rename_contigs(contig_path)