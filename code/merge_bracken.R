library(tidyverse)
source(here::here("code/util_funcs.R"))

args = commandArgs(trailingOnly=TRUE)

# Parse the GTDB database inspection
gtdb_kraken_tax <- parse_kraken_inspect(args[1])

# Further cleanup the inspection and make it match refseq
gtdb_kraken_taxa_clean <- gtdb_kraken_tax %>% 
  dplyr::select(taxID,D:S) %>% 
  mutate(across(D:S, ~str_remove(.x,"^[a-z]_")),
         taxID = paste0("gtdb_",taxID)) %>% 
  unite("taxonomy",D:S,sep = "; ",remove = FALSE) %>% 
  mutate(taxonomy = paste0("r__cellular organisms; ", taxonomy)) %>% 
  dplyr::rename(Domain = "D"
                , Phylum = "P"
                , Class = "C"
                , Order = "O"
                , Family = "F"
                , Genus = "G"
                , Species = "S"
                , taxonomy_id = "taxID"
  ) %>% 
  write_tsv("data/reference/gtdb_kraken_tax.tsv")

print("Inspected gtdb Kraken database.")

# Parse the refseq database inspection file
refseq_kraken_tax <- parse_kraken_inspect(args[2])


refseq_kraken_taxa_clean <- refseq_kraken_tax %>% 
  dplyr::select(taxID, lineage,R:ncol(.)) %>% 
  mutate(across(R:ncol(.), ~str_remove(.x,"^[a-z]_")),
         taxID = paste0("refseq_",taxID)) %>% 
  mutate(R1 = str_replace(R1,pattern = "^r1_",replacement = "r__"),
         D = paste0("d__",D),
         P = paste0("p__",P),
         C = paste0("c__",C),
         O = paste0("o__",O),
         F = paste0("f__",F),
         G = paste0("g__",G),
         S = paste0("s__",S)
  ) %>% 
  dplyr::rename(Root = "R1"
                , Domain = "D"
                , Phylum = "P"
                , Class = "C"
                , Order = "O"
                , Family = "F"
                , Genus = "G"
                , Species = "S"
                , full_lineage = "lineage"
                , taxonomy_id = "taxID"
  ) %>% 
  unite("taxonomy",c("Root","Domain","Phylum","Class","Order","Family","Genus","Species"),sep = "; ",remove = FALSE,na.rm = TRUE) %>% 
  write_tsv("data/reference/refseq_kraken_tax.tsv")

print("Inspected refseq Kraken database.")


#Find files to combine and make table of them
gtdb_files <-  data.frame(file = list.files("data/kraken2_uniq_gtdb",full.names = TRUE)) %>%
  filter(!str_detect(file,"_for_bracken.txt"),
         str_detect(file,"_bracken.txt")) %>%
  mutate(sample = str_remove(file, "_bracken.txt"),
         sample = str_remove(sample, ".*data/kraken2_uniq_gtdb/"),
         source = row_number())


gtdb_bracken <- gtdb_files$file %>% map_dfr(read_tsv, .id = "source",num_threads = as.numeric(args[5]), show_col_types = FALSE) %>% 
  left_join(gtdb_files %>% dplyr::select(source, sample) %>% mutate(source = as.character(source))) %>% 
  dplyr::select(sample, everything(), -source)

print("Imported gtdb bracken files.")

wide_gtdb_bracken <- gtdb_bracken %>% 
  #filter(new_est_reads > 2) %>% 
  dplyr::select(taxonomy_id,sample,new_est_reads) %>% 
  mutate(taxonomy_id = paste0("gtdb_",taxonomy_id)) %>% 
  pivot_wider(names_from = "sample", values_from = new_est_reads,values_fill = 0) %>% 
  left_join(gtdb_kraken_taxa_clean %>% dplyr::select(taxonomy_id,taxonomy)) %>% 
  relocate(taxonomy_id,taxonomy)

print("Processed gtdb bracken files.")

#Find files to combine and make table of them
refseq_files <- data.frame(file = list.files("data/kraken2_uniq_refseq",full.names = TRUE)) %>%
  filter(!str_detect(file,"_for_bracken.txt"),
         str_detect(file,"_bracken.txt")) %>%
  mutate(sample = str_remove(file, "_bracken.txt"),
         sample = str_remove(sample, ".*data/kraken2_uniq_refseq/"),
         source = row_number())


refseq_bracken <- refseq_files$file %>% map_dfr(read_tsv, .id = "source", num_threads = as.numeric(args[5]), show_col_types = FALSE) %>% 
  left_join(refseq_files %>% dplyr::select(source, sample) %>% mutate(source = as.character(source))) %>% 
  dplyr::select(sample, everything(), -source)

print("Imported refseq bracken files.")

wide_refseq_bracken <- refseq_bracken %>% 
  #filter(new_est_reads > 2) %>%
  dplyr::select(taxonomy_id,sample,new_est_reads) %>% 
  mutate(taxonomy_id = paste0("refseq_",taxonomy_id)) %>% 
  pivot_wider(names_from = "sample", values_from = new_est_reads,values_fill = 0) %>% 
  left_join(refseq_kraken_taxa_clean %>% dplyr::select(taxonomy_id,taxonomy)) %>% 
  relocate(taxonomy_id,taxonomy) %>% 
  filter(!str_detect(taxonomy, "d__Bacteria"),
         !str_detect(taxonomy, "d__Archaea")) 

print("Processed refseq bracken files.")

combined_braken <- bind_rows(wide_gtdb_bracken,wide_refseq_bracken) %>% 
  write_tsv(args[3])
