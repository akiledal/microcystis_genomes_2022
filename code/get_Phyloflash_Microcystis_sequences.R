# Set working directory to project root directory
setwd(here::here())

library(tidyverse)

# Read in phyloFlash classifications of the assembled 16S sequences
spades_classifications <- list.files(path = "data/phyloFlash", pattern = "phyloFlash.extractedSSUclassifications.csv", full.names = T, recursive = T) %>% 
  map_df(read_csv)

# Read in sequencingID to strainID map file
seq_to_strain <- read_tsv("data/strain_to_seq_sample_map.tsv")

# Function for reading in fasta to dataframe
fasta_df <- function(path){
  Biostrings::readDNAStringSet(path) %>% 
    data.frame(seq = .) %>% 
    rownames_to_column("header")
}

# Find and read in all classification files
spades_16S_seqs <- list.files(path = "data/phyloFlash", pattern = ".all.final.fasta", full.names = T, recursive = T) %>% 
  map_df(fasta_df) %>% 
  mutate(OTU = str_remove(header,"_[0-9,\\.]*$")) %>% 
  left_join(spades_classifications)

# Get only the Microcystis sequences
microcystis_seqs <- spades_16S_seqs %>% 
  filter(str_detect(taxonomy,"Microcystis")) %>% 
  select(header,seq) %>% 
  mutate(sequencingID = str_remove(header, "\\.PFspades.*")) %>% 
  left_join(seq_to_strain) %>%
  mutate(new_header = glue::glue("{strainID}__{sequencingID}")) %>% 
  select(header = "new_header", seq, sample="sequencingID")

microcystis_biostrings <- Biostrings::DNAStringSet(microcystis_seqs$seq)
names(microcystis_biostrings) <- microcystis_seqs$header

# Export the Microcystis sequences
Biostrings::writeXStringSet(microcystis_biostrings,"data/phyloFlash/microcystis_seqs.fasta")


# Export individual microcystis seqs

# microcystis_seqs_w_sample <- microcystis_seqs %>% 
#   mutate(sample = str_remove(header, "\\.PFspades.*"))

export_indiv_seq <- function(sample_str){
  export_seq <- microcystis_seqs %>% 
    filter(sample == sample_str)
  
  indiv_seq <- Biostrings::DNAStringSet(export_seq$seq)
  names(indiv_seq) <- export_seq$header
  
  # Export the Microcystis sequences
  Biostrings::writeXStringSet(indiv_seq, glue::glue("data/phyloFlash/{sample_str}/{sample_str}_microcystis_16S.fasta"))
  Biostrings::writeXStringSet(indiv_seq, glue::glue("data/phyloFlash/read_mapping/{sample_str}_microcystis_16S.fasta"))
}

sapply(microcystis_seqs$sample,export_indiv_seq)
