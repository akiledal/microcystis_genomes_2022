setwd(here::here())
library(tidyverse)

#####
# Tidy the marker gene variants determined by variant calling from mapped reads
#####

microcystis_marker_genes <- list.files("data/marker_genes", "_consensus.fasta",full.names = TRUE, recursive = T) %>% 
  data.frame(path = .) %>% 
  mutate(gene = str_remove(path,".*/") %>% str_remove_all("_.*"),
         sample = str_remove(path,"_consensus.fasta") %>% str_remove(glue::glue(".*/{gene}_")),
         seq = Biostrings::readDNAStringSet(path) %>% as.character())

genes <- microcystis_marker_genes %>% 
  select(gene) %>% 
  distinct() %>% 
  pull("gene")

write_fastas <- function(gene_name){
  to_export <- microcystis_marker_genes %>% 
    filter(gene == gene_name) %>% 
    mutate(sample = paste0(gene,"_",sample,"_consensus")) %>% 
    select(sample, seq)
  
  microcystis_biostrings <- Biostrings::DNAStringSet(to_export$seq)
  names(microcystis_biostrings) <- to_export$sample
  
  dir.create("data/marker_genes/identified_variants")
  Biostrings::writeXStringSet(microcystis_biostrings,glue::glue("data/marker_genes/identified_variants/{gene_name}_called_variants.fasta"))
  
}

sapply(genes,write_fastas)


#####
# Tidy the marker gene variants determined by assembling marker genes with reads mapping to reference sequences
#####

assembled_microcystis_marker_genes <- list.files("data/marker_genes", "_assembly.fasta",full.names = TRUE, recursive = T) %>% 
  data.frame(path = .)

fasta_df <- function(path){
  Biostrings::readDNAStringSet(path) %>% 
    data.frame(seq = .) %>% 
    mutate(gene = str_remove(path,".*/") %>% str_remove_all("_assembly.fasta"),
           sample = str_remove(path, "data/marker_genes/") %>% str_remove("/.*")) %>% 
    rownames_to_column("header")
}

assembled_markers <- map_df(assembled_microcystis_marker_genes$path,fasta_df) %>% 
  group_by(gene, sample) %>% 
  mutate(row_number = row_number(),
         new_header = glue::glue("{gene}_{sample}_assembled_{row_number}")) %>% 
  ungroup()

genes <- assembled_markers %>% 
  select(gene) %>% 
  distinct() %>% 
  pull("gene")

write_fastas <- function(gene_name){
  to_export <- assembled_markers %>% 
    filter(gene == gene_name) %>% 
    select(new_header, seq)
  
  microcystis_biostrings <- Biostrings::DNAStringSet(to_export$seq)
  names(microcystis_biostrings) <- to_export$new_header
  
  dir.create("data/marker_genes/identified_variants")
  Biostrings::writeXStringSet(microcystis_biostrings,glue::glue("data/marker_genes/identified_variants/{gene_name}_assembled.fasta"))
  
}

sapply(genes,write_fastas)