'Prepare ghostKOALA input.

Usage:
  prepare_for_ghostKOALA.R [options]

Options:
  -i DIR --in_dir=DIR         Directory with prodigal genes
  -t TAX --mag-tax=TAX        GTDBtk MAG taxonomy file
  -d DREP --drep_bins=DREP    Representative bins determined by drep
  -m MAP --mag-map=MAP        Output file mapping bin name to ID for ghostKOALA 
  -o FAA --combined-fasta     Combined fasta for ghostKOALA, with headers formatted for KEGGdecoder
  -h --help                   Show this help screen.
' -> doc

library(docopt)

# To process actual command line input
arguments <- docopt(doc)
#print(arguments)

# # for testing interactively
# arguments <- docopt(doc, args = c("-i data/prodigal_mag -t data/gtdb/gtdbtk.bac120.summary.tsv -d data/drep/data_tables/Wdb.csv -m data/kegg/mag_name_map.tsv -o data/kegg/mag_prodigal_combined.faa"))
# print(arguments)


library(tidyverse)

setwd(here::here()) #set working directory to project root
# List prodigal gene calls for each bin
amino_acid_fastas <- list.files(arguments$in_dir,pattern = "*.faa",full.names = TRUE)

# Import bin taxonomy
mag_tax <- read_tsv(arguments$mag_tax) %>%
  dplyr::select(-red_value)

# Import drep results
drep_bins <- read_csv(arguments$drep_bins) %>% pull("genome") %>% str_remove(".fasta")

# Make map file from old bin names to simplified names matching KEGGdecoder expectations
mag_map <- data.frame(mag_fp = amino_acid_fastas) %>%
  mutate(old_mag_name = str_remove(mag_fp,paste0(arguments$in_dir,"/")) %>% str_remove(".faa")) %>%
  arrange(old_mag_name) %>%
  mutate(kegg_decoder_name = glue::glue("MAG{row_number()}")) %>% 
  left_join(mag_tax %>% dplyr::select(user_genome, classification),by = c("old_mag_name" = "user_genome")) %>%
  write_tsv(arguments$mag_map)

# Function for reading in .faa 
mod_fasta <- function(fasta_path){

  aa <- data.frame(seq = Biostrings::readAAStringSet(fasta_path)) %>%
    rownames_to_column("header") %>% mutate(mag_fp = fasta_path)
}

# Replace header with simplified header for KEGGdecoder
mod_fastas <- lapply(amino_acid_fastas,mod_fasta) %>%
  bind_rows() %>%
  left_join(mag_map) %>%
  group_by(kegg_decoder_name) %>%
  mutate(header = glue::glue("{kegg_decoder_name}_{row_number()}")) %>%
  ungroup() %>%
  dplyr::select(header,seq)

fasta_export <- mod_fastas$seq
names(fasta_export) <-mod_fastas$header

print("Headers cleaned, writing concatenated .faa")

Biostrings::AAStringSet(fasta_export,use.names = TRUE) %>% 
  Biostrings::writeXStringSet(arguments$combined_fasta)