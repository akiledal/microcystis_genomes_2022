'Prepare ghostKOALA input.

Usage:
  prepare_for_ghostKOALA.R [options]

Options:
  -i DIR --in_dir=DIR         Directory with prodigal gene calls
  -o FAA --combined-fasta     Combined fasta for ghostKOALA, with headers formatted for KEGGdecoder
  -h --help                   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
# arguments <- docopt(doc, args = c("--prodigal_dir=data/prodigal_mag"))
# print(arguments)


library(tidyverse)

setwd(here::here()) #set working directory to project root
# List prodigal gene calls for each bin

amino_acid_fastas <- list.files(arguments$in_dir,pattern = "*.faa",full.names = TRUE)

mag_map <- data.frame(mag_fp = amino_acid_fastas) %>% 
  mutate(assembly = str_remove(mag_fp,paste0(arguments$in_dir,"/")) %>% str_remove(".faa"))

mod_fasta <- function(fasta_path){
  
  aa <- data.frame(seq = Biostrings::readAAStringSet(fasta_path)) %>%
    rownames_to_column("header") %>% mutate(mag_fp = fasta_path)
}

mod_fastas <- lapply(amino_acid_fastas,mod_fasta) %>%
  bind_rows() %>%
  left_join(mag_map) %>%
  group_by(assembly) %>%
  mutate(header = glue::glue("{assembly}_{row_number()}")) %>%
  ungroup() %>%
  dplyr::select(header,seq)

fasta_export <- mod_fastas$seq
names(fasta_export) <-mod_fastas$header

Biostrings::AAStringSet(fasta_export,use.names = TRUE) %>% Biostrings::writeXStringSet(arguments$combined_fasta)