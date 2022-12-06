# This file contains useful utility scripts used in the metagenome analysis


# This function processes a kraken2 database inspection and produces a
#   dataframe of taxonomy information in a format more easily used with
#   tidyverse tools. It relies heavily on Pavian, from the same from that produced
#   Kraken2. A similar approach could be used for each kraken2 report, but it is much
#   slower than just doing this once for all possible taxa.
parse_kraken_inspect <- function(inspect_file){

  # Define column names, since no headers in the inspection file
  cols <- c("ra_and_subs","reads_and_subs","reads","rank","id","tax")
  
  # Read in inspection file and produce table with tax levels in the proper order
  tax_levels <- read_tsv(inspect_file, col_names = cols,comment = "#") %>% 
    pull("rank") %>% 
    unique() %>% 
    data.frame(level = .) %>% 
    mutate(toplevel = str_extract(level,"[A-Z]"),
           sublevel = str_extract(level,"[0-9]"),
           sublevel = if_else(is.na(sublevel), 0, as.numeric(sublevel)))
  
  # Order the taxonomy levels for pavian's importer
  tax_lev_df <- data.frame(toplevel = c("R", "D","K", "P", "C", "O", "F", "G", "S")) %>% 
    left_join(tax_levels) %>% filter(!is.na(level))
  
  #Read in the report with all taxonomic levels
  full_report <- pavian::read_report2(inspect_file,
                                      collapse = FALSE,
                                      add_taxRank_columns = TRUE,
                                      keep_taxRanks = tax_lev_df$level)
  
  # Cleanup the report, only reports species level annotations, but this could be 
  #   easily modified by changing the filtering below.
  full_clean_report <- full_report %>% 
    filter(taxRank == "S") %>% 
    mutate(name = str_remove(name,"s_")
    ) %>% 
    dplyr::select(n_read = "taxonReads",
                  taxID,
                  name,
                  lineage = "taxLineage",
                  any_of(tax_lev_df$level))
  
  return(full_clean_report)
}


# Plot abundance of a particular sample in all concrete metagenome samples and negative controls

plot_metag_abund <- function(tax_str,summarize, transform = "clr",decontam = FALSE){
  
  metadata <- read_tsv("data/metadata.tsv") %>% 
    mutate(sample_wo_replicate = str_remove(sample, "[a-z]"),
           replicate = str_remove(sample, "[A-Z,0-9]*")) %>% 
    column_to_rownames("sample") %>% 
    filter(type == "concrete" | neg_control == TRUE,
           !str_detect(rownames(.), "NEG")) %>% 
    mutate(ASR_logical = if_else(neg_control == TRUE, "neg_control", as.character(asr_logical)))
  
  
table <- if_else(decontam,"results/decontaminated_bracken_count_table.tsv","results/braken_count_table.tsv")
  
otu_table <- read_tsv(table) %>% 
    #filter(!taxonomy_id %in% contams$taxonomy_id) %>%
    dplyr::select(-taxonomy) %>% 
    column_to_rownames("taxonomy_id") %>% 
    .[,rownames(metadata)]
  
  taxonomy <- read_tsv("results/braken_count_table.tsv") %>% 
    dplyr::select(taxonomy, taxonomy_id)
  
  
  otu_table %>% 
    microbiome::transform(transform) %>% as.data.frame() %>% 
    rownames_to_column("taxonomy_id") %>% 
    pivot_longer(-taxonomy_id, names_to = "sample", values_to = "abund") %>% 
    left_join(taxonomy) %>% 
    filter(str_detect(taxonomy,tax_str)) %>%
    group_by(sample) %>% summarise(abund = sum(abund)) %>% 
    left_join(metadata %>% rownames_to_column("sample")) %>% 
    ggplot(aes(sample_wo_replicate,abund, color = location, notes =ident_info)) + 
    geom_point() + 
    facet_grid(~ ASR_logical,scales = "free_x") + 
    labs(x = NULL) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0))
}
