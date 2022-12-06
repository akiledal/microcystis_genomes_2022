#!/usr/bin/env Rscript

library(tidyverse)

files <-  data.frame(file = list.files("data/kraken2/",full.names = TRUE)) %>%
  filter(str_detect(file,"_brackenReport.txt")) %>%
  mutate(sample = str_remove(file, "_brackenReport.txt"),
         sample = str_remove(sample, "data/kraken2/"))



combined <- data.frame()
cols <- c("ra_and_subs","reads_and_subs","reads","rank","NCBI_id","tax")


for (sample_id in files$sample){

########
##Read in GTDB classified report
#######

kraken <- read_tsv(paste0("data/kraken2/",sample_id,"_brackenReport.txt"),
                   col_names = cols,trim_ws = F)

ranks_and_num <- unique(kraken$rank) %>%
  data.frame(rank = .) %>%
  mutate(num_rank = row_number())

tax_split <- ranks_and_num %>% filter(rank != "U") %>% pull(rank)

#Fill in higher level taxonomy and calculate relative abundances without sub-taxa
kraken_gtdb <- kraken %>% separate(tax,into = tax_split,sep = "  ") %>%
  mutate(across(everything(), ~na_if(.x,""))) %>%
  left_join(ranks_and_num) %>%
  fill(any_of(tax_split)) %>%
  mutate(taxa = row_number()) %>%
  pivot_longer(tax_split,names_to = "row_rank",values_to = "row_tax") %>%
  left_join(ranks_and_num %>% rename(row_rank = "rank",row_num_rank = "num_rank")) %>%
  filter(row_num_rank <= num_rank) %>%
  pivot_wider(id_cols = c("taxa","ra_and_subs","reads_and_subs","reads","rank","NCBI_id"), names_from = "row_rank",values_from = "row_tax") %>%
  mutate(rel.abund = reads/sum(reads)) %>%
  mutate(R = paste0("r__",R),
         taxonomy = paste(R,D,P,C,O,F,G,S,sep = "|"),
         taxonomy = str_remove(taxonomy,"[NA|]+$"),
         sample = sample_id) %>%
  select(sample,reads,taxonomy,R,D,P,C,O,F,G,S,rank, NCBI_id)

combined <- combined %>% bind_rows(kraken_gtdb)

}

gtdb_wide <- combined %>% pivot_wider(id_cols = c(taxonomy,R,D,P,C,O,F,G,S,rank, NCBI_id),names_from = sample, values_from = reads,values_fill = 0)


########
##Read in refseq classified report
#######

combined <- data.frame()

for (sample_id in files$sample){

kraken <- read_tsv(paste0("data/kraken2_euk/",sample_id,"_brackenReport.txt"), col_names = cols,trim_ws = F)

ranks_and_num <- unique(kraken$rank) %>%
  data.frame(rank = .) %>%
  mutate(num_rank = row_number())

tax_split <- ranks_and_num %>% filter(rank != "U") %>% pull(rank)

#Fill in higher level taxonomy and calculate relative abundances without sub-taxa
kraken_refseq <- kraken %>% separate(tax,into = tax_split,sep = "  ") %>%
  mutate(across(everything(), ~na_if(.x,""))) %>%
  left_join(ranks_and_num) %>%
  fill(any_of(tax_split)) %>%
  mutate(taxa = row_number()) %>%
  pivot_longer(tax_split,names_to = "row_rank",values_to = "row_tax") %>%
  left_join(ranks_and_num %>% rename(row_rank = "rank",row_num_rank = "num_rank")) %>%
  filter(row_num_rank <= num_rank) %>%
  pivot_wider(id_cols = c("taxa","ra_and_subs","reads_and_subs","reads","rank","NCBI_id"), names_from = "row_rank",values_from = "row_tax") %>%
  mutate(rel.abund = reads/sum(reads),
         across(c(R,D,P,C,O,F,G,S),~paste(tolower(cur_column()),.x,sep = "__")),
         taxonomy = paste(R,D,P,C,O,F,G,S,sep = "|"),
         taxonomy = str_remove(taxonomy,"[NA|]+$"),
         sample = sample_id) %>%
  filter(D != "d__Bacteria") %>%
  select(sample,reads,taxonomy,R,D,P,C,O,F,G,S,rank, NCBI_id)

combined <- combined %>% bind_rows(kraken_refseq)

}

refseq_wide <- combined %>% pivot_wider(id_cols = c(taxonomy,R,D,P,C,O,F,G,S,rank, NCBI_id),names_from = sample, values_from = reads,values_fill = 0)

abudance <- bind_rows(gtdb_wide, refseq_wide)

#Combine bacteria and archaea from GTDB with other life from refseq


rel.abund <- abudance %>%
  mutate(across(-c(taxonomy,R,D,P,C,O,F,G,S,rank, NCBI_id), ~(.x/sum(.x)))) %>% write_tsv("results/combined_bracken.tsv")
