import os
import re
from glob import glob
from snakemake.utils import R

# Set config file path
configfile: "config.yaml"

# Where to save Snakemake report
report: "code/report/workflow.rst"

# Determine sample names from read files
samples = glob_wildcards("data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq").sample

# Target rule for producing all output
rule all:
    input: 
        expand("data/kraken2_uniq_{database}/{sample}_out.txt", database = ["refseq","gtdb_r202"], sample = samples),
        expand("data/kraken2_uniq_{database}/{sample}_brackenMpa.txt", database = ["refseq","gtdb_r202"], sample = samples),
        expand("results/metacodeR/{sample}_kraken.pdf", sample=samples),
        "data/drep/dereplicated_genomes",
        "data/gtdb",
        "results/bin_coverage.tsv",
        expand("data/prodigal_mag/{MAG_name}.gbk", MAG_name=glob_wildcards('data/combined_bins/{MAG}.fasta').MAG),
        expand("data/prodigal_assembly/{sample}.gbk", sample = samples),
        expand("data/plass/{sample}_assembly.fa",sample = samples)


# Check reads
rule multiqc_PREqc:
    input:
        f_seqs = expand("data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",sample = samples),
        r_seqs = expand("data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq", sample = samples)
    output: 
        directory("results/multiqc_PREqc")
    conda: "code/multiqc_env.yaml"
    shell:
        """
        multiqc /omics/HABs/Microcystis/Microcystis_StrainDB_ROS_Enzyme_Search/WLE_strain_isolates/qcd_reads/*/pre_fastq -o {output}

        """

rule multiqc_POSTqc:
    input:
        f_seqs = expand("data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",sample = samples),
        r_seqs = expand("data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq", sample = samples)
    output: 
        directory("results/multiqc_POSTqc")
    conda: "code/multiqc_env.yaml"
    shell:
        """
        multiqc /omics/HABs/Microcystis/Microcystis_StrainDB_ROS_Enzyme_Search/WLE_strain_isolates/qcd_reads/*/post_fastq -o {output}

        """


# Target rules for running kraken
rule calc_kraken_uniq:
    input: expand("data/kraken2_uniq_{database}/{sample}_out.txt", database = ["refseq","gtdb_r202"], sample = samples),
        expand("data/kraken2_uniq_{database}/{sample}_brackenMpa.txt", database = ["refseq","gtdb_r202"], sample = samples)

rule calc_kraken_uniq_gtdb:
    input: expand("data/kraken2_uniq_{database}/{sample}_out.txt", database = ["gtdb_r202"], sample = samples)

# Inspect Kraken2 databases, nescessary to produce nicely formatted taxonomy files for downstream analysis
rule kraken_inspect:
    input:
        db = ancient("/geomicro/data2/kiledal/references/kraken_databases/{database}")
    output: 
        inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/{database}/inspect.txt"
    conda: "code/main_env.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        kraken2-inspect --db {input.db} > {output.inspect_file}
        """

rule add_lineage_to_inspect_gtdb:
    input:
        db = ancient("/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202"),
        inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect.txt"
    output: 
        inspect_w_lineage = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect_w_lineage.txt"
    conda: "code/main_env.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        taxonkit lineage \
            {input.inspect_file} \
            --data-dir {input.db}/taxonomy \
            -i 5 \
            -o {output.inspect_w_lineage}
        """

rule add_lineage_to_inspect_refseq:
    input:
        db = ancient("/geomicro/data2/kiledal/references/kraken_databases/refseq"),
        inspect_file = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect.txt"
    output: 
        inspect_w_lineage_unformatted = temp("/geomicro/data2/kiledal/references/kraken_databases/refseq/unformatted_inspect_w_lineage.txt"),
        inspect_w_lineage = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect_w_lineage.txt"
    conda: "code/main_env.yaml"
    resources: cpus=1, mem_mb=250000, time_min=5440, mem_gb = 250
    shell:
        """
        taxonkit lineage \
            {input.inspect_file} \
            --data-dir {input.db}/taxonomy \
            -i 5 \
            -o {output.inspect_w_lineage_unformatted}

        taxonkit reformat \
            {output.inspect_w_lineage_unformatted} \
            -i 7 \
            -P \
            -o {output.inspect_w_lineage}
        """


rule kraken_database_tax_merge:
    input:
        script = "code/merge_kraken_tax.R",
        gtdb_tax_info = rules.add_lineage_to_inspect_refseq.output.inspect_w_lineage,
        refseq_tax_info = rules.add_lineage_to_inspect_gtdb.output.inspect_w_lineage
    output:
        combined_tax_info = "data/reference/kraken_tax_info_merged.tsv"
    resources: cpus=1, mem_mb=4000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        ./{input.script} -g {input.gtdb_tax_info} -r {input.refseq_tax_info} -o {output.combined_tax_info}
        """


# Run kraken2 with KrakenUniq like functionality to screen for higher confidence results
# First runs with a GTDB database to annotate bacterial and archaeal reads
rule kraken2_gtdb_w_uniq:
    input:
        f_seq = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_seq = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        db = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202"
    output:
        report = "data/kraken2_uniq_gtdb_r202/{sample}_report.txt",
        out = temp("data/kraken2_uniq_gtdb_r202/{sample}_out.txt"),
        bracken = "data/kraken2_uniq_gtdb_r202/{sample}_bracken.txt",
        bracken_report = "data/kraken2_uniq_gtdb_r202/{sample}_brackenReport.txt",
        bracken_mpa = "data/kraken2_uniq_gtdb_r202/{sample}_brackenMpa.txt",
        unclass_f = temp("data/kraken2_uniq_gtdb_r202/{sample}_unclassified_1.fasta"),
        unlcass_r = temp("data/kraken2_uniq_gtdb_r202/{sample}_unclassified_2.fasta")
    params:
        uniq_minimizer_threshold = 150
    conda: "code/main_env.yaml"
    resources: cpus=8, mem_mb=252000, time_min=1440, mem_gb = 250
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            --db {input.db} \
            --minimum-hit-groups 3 \
            --unclassified-out data/kraken2_uniq_gtdb_r202/{wildcards.sample}_unclassified#.fasta \
            --paired {input.f_seq} {input.r_seq}

        echo "Kraken complete, filtering..."

        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > data/kraken2_uniq_gtdb_r202/{wildcards.sample}_for_bracken.txt

        echo "Filtering complete, running bracken..."

        bracken -d {input.db} -i data/kraken2_uniq_gtdb_r202/{wildcards.sample}_for_bracken.txt -o {output.bracken} -w {output.bracken_report}

        ./code/kreport2mpa.py -r {output.bracken_report} -o {output.bracken_mpa} --percentages

        echo "Bracken complete. Quitting."
        """

# Any reads not annotated with the GTDB database are then annotated with a RefSeq database
rule kraken2_refseq_w_uniq: ##Run kraken2
    input:
        f_seq = "data/kraken2_uniq_gtdb_r202/{sample}_unclassified_1.fasta",
        r_seq = "data/kraken2_uniq_gtdb_r202/{sample}_unclassified_2.fasta",
        db = "/geomicro/data2/kiledal/references/kraken_databases/refseq"
    output:
        report = "data/kraken2_uniq_refseq/{sample}_report.txt",
        out = temp("data/kraken2_uniq_refseq/{sample}_out.txt"),
        bracken = "data/kraken2_uniq_refseq/{sample}_bracken.txt",
        bracken_report = "data/kraken2_uniq_refseq/{sample}_brackenReport.txt",
        bracken_mpa = "data/kraken2_uniq_refseq/{sample}_brackenMpa.txt"
    params:
        uniq_minimizer_threshold = 150
    conda: "code/main_env.yaml"
    resources: cpus=8, mem_mb=250000, time_min=1440, mem_gb = 250
    shell:
        """
        kraken2 \
            --threads {resources.cpus} \
            --report-minimizer-data \
            --report {output.report} \
            --output {output.out} \
            --db {input.db} \
            --minimum-hit-groups 3 \
            --paired {input.f_seq} {input.r_seq}

        awk '{{ if($5 >= {params.uniq_minimizer_threshold}) {{ print }}}}' {output.report} | cut --complement -f4,5 > data/kraken2_uniq_refseq/{wildcards.sample}_for_bracken.txt

        bracken -d {input.db} -i data/kraken2_uniq_refseq/{wildcards.sample}_for_bracken.txt -o {output.bracken} -w {output.bracken_report}

        ./code/kreport2mpa.py -r {output.bracken_report} -o {output.bracken_mpa} --percentages
        """

# Combine the annotations and produce count table
rule kraken_summarize:
    input:
        kraken_results = expand("data/kraken2_uniq_{database}/{sample}_bracken.txt", database = ["refseq","gtdb_r202"], sample = samples),
        gtdb_inspect = "/geomicro/data2/kiledal/references/kraken_databases/gtdb_r202/inspect.txt",
        refseq_inspect = "/geomicro/data2/kiledal/references/kraken_databases/refseq/inspect.txt"
    output:
        table = "results/braken_count_table.tsv",
        gtdb_tax = "data/reference/gtdb_kraken_tax.tsv",
        refseq_tax = "data/reference/refseq_kraken_tax.tsv"
    resources: cpus=4, mem_mb=100000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
         which R
         R --vanilla << 'RSCRIPT'

        library(tidyverse)

        source(here::here("code/util_funcs.R"))
        setwd(here::here())

        # Parse the GTDB database inspection
        gtdb_kraken_tax <- parse_kraken_inspect("{input.gtdb_inspect}")

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
        write_tsv("{output.gtdb_tax}")

        # Parse the refseq database inspection file
        refseq_kraken_tax <- parse_kraken_inspect("{input.refseq_inspect}")

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
        write_tsv("{output.refseq_tax}")


        #Find files to combine and make table of them
        gtdb_files <-  data.frame(file = list.files("data/kraken2_uniq_gtdb_r202",full.names = TRUE)) %>%
        filter(!str_detect(file,"_for_bracken.txt"),
                str_detect(file,"_bracken.txt")) %>%
        mutate(sample = str_remove(file, "_bracken.txt"),
                sample = str_remove(sample, ".*data/kraken2_uniq_gtdb_r202/"),
                source = row_number())


        gtdb_bracken <- gtdb_files$file %>% map_dfr(read_tsv, .id = "source",num_threads = 4, show_col_types = FALSE) %>% 
            left_join(gtdb_files %>% dplyr::select(source, sample) %>% mutate(source = as.character(source))) %>% 
            dplyr::select(sample, everything(), -source)


        wide_gtdb_bracken <- gtdb_bracken %>% 
            filter(new_est_reads > 2) %>% 
            dplyr::select(taxonomy_id,sample,new_est_reads) %>% 
            mutate(taxonomy_id = paste0("gtdb_",taxonomy_id)) %>% 
            pivot_wider(names_from = "sample", values_from = new_est_reads,values_fill = 0) %>% 
            left_join(gtdb_kraken_taxa_clean %>% dplyr::select(taxonomy_id,taxonomy)) %>% 
            relocate(taxonomy_id,taxonomy)


        #Find files to combine and make table of them
        refseq_files <- data.frame(file = list.files("data/kraken2_uniq_refseq",full.names = TRUE)) %>%
            filter(!str_detect(file,"_for_bracken.txt"),
                    str_detect(file,"_bracken.txt")) %>%
            mutate(sample = str_remove(file, "_bracken.txt"),
                    sample = str_remove(sample, ".*data/kraken2_uniq_refseq/"),
                    source = row_number())


        refseq_bracken <- refseq_files$file %>% map_dfr(read_tsv, .id = "source", num_threads = 4, show_col_types = FALSE) %>% 
            left_join(refseq_files %>% dplyr::select(source, sample) %>% mutate(source = as.character(source))) %>% 
            dplyr::select(sample, everything(), -source)


        wide_refseq_bracken <- refseq_bracken %>% 
            filter(new_est_reads > 2) %>%
            dplyr::select(taxonomy_id,sample,new_est_reads) %>% 
            mutate(taxonomy_id = paste0("refseq_",taxonomy_id)) %>% 
            pivot_wider(names_from = "sample", values_from = new_est_reads,values_fill = 0) %>% 
            left_join(refseq_kraken_taxa_clean %>% dplyr::select(taxonomy_id,taxonomy)) %>% 
            relocate(taxonomy_id,taxonomy) %>% 
            filter(!str_detect(taxonomy, "d__Bacteria"),
                    !str_detect(taxonomy, "d__Archaea")) 

        combined_braken <- bind_rows(wide_gtdb_bracken,wide_refseq_bracken) %>% 
            write_tsv("{output.table}")

        'RSCRIPT'
        """

# Target rule to make all the metacodeR plots
rule plot_metacoders:
    input: expand("results/metacodeR/{sample}_kraken.pdf", sample=samples)

# Produce metacodeR figures for each sample summarizing the community compostion and relative abundance
rule metacodeR:
    input:
        combined = "results/braken_count_table.tsv"
        #metadata = "data/metadata.tsv"
    output: "results/metacodeR/{sample}_kraken.pdf"
    resources: cpus=1, mem_mb=8000, time_min=60
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        mkdir -p results/metacodeR

        R --vanilla << 'RSCRIPT'

        library(metacoder)
        library(tidyverse)

        rel.abund <- read_tsv("{input.combined}") %>% 
            mutate(across(-c("taxonomy_id", "taxonomy"), ~ .x / sum(.x)))

        sample_abund <- rel.abund %>%
        dplyr::select(`{wildcards.sample}`,taxonomy) %>%
        filter(`{wildcards.sample}` > 0)


        obj <- parse_tax_data(sample_abund, class_cols = "taxonomy", class_sep = "; ",
                            class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                            class_regex = "^(.+)__(.*)$")

        #converting to relative abundance
        obj$data$tax_data <- calc_obs_props(obj, "tax_data")

        #summing per-taxon counts
        obj$data$tax_abund <- calc_taxon_abund(obj, "tax_data")

        obj$data$tax_abund <- obj$data$tax_abund %>% rename(rel_abund = "{wildcards.sample}")

        #obj2 <- obj %>% filter_taxa(taxon_ranks == "g", supertaxa = TRUE, reassign_obs = TRUE, drop_obs = TRUE)

        (tree <- obj %>%
            metacoder::filter_taxa(taxon_ranks == "g", supertaxa = TRUE,reassign_obs = TRUE, subtaxa = FALSE) %>% # subset to the genus rank
            metacoder::filter_obs("tax_abund", rel_abund > 0.005,drop_taxa = TRUE,supertaxa = TRUE,reassign_obs = TRUE) %>% #filter to only plot taxa that represent at least 0.5% of the community
            heat_tree(node_label = taxon_names,
                    node_size = rel_abund,
                    node_size_range = c(0.00175,0.045),
                    node_label_size_range = c(.015,.025),
                    node_size_axis_label = "OTU count",
                    initial_layout = "reingold-tilford", layout = "davidson-harel",
                    overlap_avoidance = 10,
                    node_label_max = 75,
                    node_color = rel_abund,
                    node_color_range = c("gray","gray","gray"),
                    node_color_axis_label = "Relative abundance")
        )
        ggsave(plot = tree,filename = "{output}", device = cairo_pdf(), width = 12, height = 12, dpi = 600)

        'RSCRIPT'
        """

# Run drep to determine unique/distinct bins, and filter low quality bins
rule drep:
    output:
        main_dir = directory("data/drep"),
        MAGs= directory("data/drep/dereplicated_genomes"),
        table = "data/drep/data_tables/Wdb.csv"
    conda: "code/drep.yaml"
    resources: cpus=8, mem_mb=250000, time_min=2880,
    shell:
        """
        dRep dereplicate {output.main_dir} -g data/combined_bins/*.fasta
        """

# Download and extract reference files for the GTDB toolkit
rule download_gtdbtk_refs:
    output:
        dir = directory("/geomicro/data2/kiledal/references/gtdbtk"),
        tar = "/geomicro/data2/kiledal/references/gtdbtk/gtdbtk_data.tar.gz",
        release202 = directory("/geomicro/data2/kiledal/references/gtdbtk/release202")
    resources: cpus=1, mem_mb=8000, time_min=2880, mem_gb = 8
    shell:
        """
        mkdir -p {output.dir}
        cd {output.dir}
        wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_data.tar.gz
        tar xzf gtdbtk_data.tar.gz
        """

# Run the GTDB toolkit on bins
rule gtdbtk:
    input:
        bins = "data/combined_bins",
        refs = "/geomicro/data2/kiledal/references/gtdbtk/release202"
    output: 
        dir = directory("data/gtdb"),
        bac_tax = "data/gtdb/gtdbtk.bac120.summary.tsv "
    conda: "code/gtdbtk.yaml"
    resources: cpus=1, mem_mb=500000, time_min=2880, mem_gb = 500
    shell:
        """
        GTDBTK_DATA_PATH={input.refs}

        gtdbtk classify_wf --extension fasta --genome_dir {input.bins} --out_dir {output.dir} --cpus {resources.cpus}
        """

# Map QCd reads to MAGs using coverm which uses the minimap2 aligner.
rule coverm:
    input:
        bins = "data/combined_bins"
    output: "results/bin_coverage.tsv"
    conda: "code/coverm_env.yaml"
    resources: cpus=16, mem_mb=250000, time_min=2880
    shell:
        """
        coverm genome \
            -t {resources.cpus} \
            -m relative_abundance mean covered_bases variance length \
            --min-covered-fraction 0 \
            --coupled data/qcd_reads/*deduped_trimmed_screened.fq \
            --genome-fasta-files {input.bins}/*.fasta \
            -o {output}
        """

rule coverm_drepd:
    input:
        bins = "data/drep/dereplicated_genomes"
    output: "results/deduped_bin_coverage.tsv"
    conda: "code/coverm_env.yaml"
    resources: cpus=16, mem_mb=250000, time_min=2880
    shell:
        """
        coverm genome \
            -t {resources.cpus} \
            -m relative_abundance mean covered_bases variance length \
            --min-covered-fraction 0 \
            --coupled data/qcd_reads/*deduped_trimmed_screened.fq \
            --genome-fasta-files {input.bins}/*.fasta \
            -o {output}
        """


# Use prodigal to call genes in assemblies
rule prodigal_assembly:
    input:
        assembly = "data/assembly/{sample}/Megahit_meta-sensitive_out/final.contigs.fa"
    output: 
        gbk="data/prodigal_assembly/{sample}.gbk",
        faa="data/prodigal_assembly/{sample}.faa"
    log:
        "logs/prodigal_assembly/{sample}.log"
    conda: "code/main_env.yaml"
    resources: cpus=1, mem_mb=10000, time_min=1440
    shell:
        """
        prodigal -i {input.assembly} -o {output.gbk} -a {output.faa} 2> {log}
        """

# Target rule for running prodigal on assemblies and bins
rule run_prodigal:
    input: 
        expand("data/prodigal_mag/{MAG_name}.gbk", MAG_name=glob_wildcards('data/combined_bins/{MAG}.fasta').MAG),
        expand("data/prodigal_assembly/{sample}.gbk", sample = samples)

# Use prodigal to call genes in MAGs.
rule prodigal_mags:
    input:
        mag = "data/combined_bins/{MAG}.fasta"
    output:
        touch="data/prodigal_mag/touches/{MAG}",
        faa="data/prodigal_mag/{MAG}.faa"
    params:
        gbk="data/prodigal_mag/{MAG}.gbk",
        faa="data/prodigal_mag/{MAG}.faa",
    log:
        "logs/prodigal_mag/{MAG}.log"
    conda: "code/main_env.yaml"
    resources: cpus=1, mem_mb=10000, time_min=1440
    shell:
        """
        touch {output.touch}
        prodigal -i {input.mag} -o {params.gbk} -a {params.faa} 2> {log} || true
        """

rule ghostKOALA_bin_prep:
    input:
        genes= expand("data/prodigal_mag/touches/{MAG_name}", MAG_name=glob_wildcards('data/combined_bins/{MAG}.fasta').MAG),
        bin_taxonomy = "data/gtdb/gtdbtk.bac120.summary.tsv",
        drep = "data/drep/data_tables/Wdb.csv"
    output:
        mag_map = "data/kegg/mag_name_map.tsv",
        merged_genes = "data/kegg/mag_prodigal_combined.faa"
    params:
        mag_prodigal_genes_dir = "data/prodigal_mag",
        script = "code/prepare_MAG_genes_for_ghostKOALA.R"
    log: "logs/prep_for_ghostKOALA.log"
    resources: cpus = 1, mem_mb = 4000, time_min = 20
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        Rscript {params.script} \
            -i {params.mag_prodigal_genes_dir} \
            -t {input.bin_taxonomy} \
            -d {input.drep} \
            -m {output.mag_map} \
            -o {output.merged_genes}
        """

rule ghostKOALA_assembly_prep:
    input:
        expand("data/prodigal_assembly/{sample}.faa",sample = samples)
    output:
        merged_genes = "data/kegg/assembly_prodigal_combined.faa"
    params:
        mag_prodigal_genes_dir = "data/prodigal_assembly",
        script = "code/prepare_assembly_genes_for_ghostKOALA.R"
    log: "logs/prep_for_ghostKOALA.log"
    resources: cpus = 1, mem_mb = 4000, time_min = 20
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        Rscript {params.script} \
            -i {params.mag_prodigal_genes_dir} \
            -o {output.merged_genes}
        """

# Target rule for running plass amino acid assemblies
rule run_plass:
    input: expand("data/plass/{sample}_assembly.fa",sample = samples)

# Run amino-acid level assembly of raw reads using PLASS. Note that plass requires a large temp directory. 
rule plass:
    input:
        #f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        #r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq"
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq"
    output: 
        assembly = "data/plass/{sample}_assembly.fa",
        temp_folder = temp("data/tmp/plass/{sample}/")
    log: 
        "logs/plass/{sample}.log"
    benchmark:
        "benchmarks/plass/{sample}.txt"
    conda: "code/plass.yaml"
    resources: cpus=12, mem_mb=40000, time_min=1440
    shell:
        """
        plass assemble \
        --threads {resources.cpus} \
        {input.f_reads} \
        {input.r_reads} \
        {output.assembly} \
        {output.temp_folder} 2> {log}
        """

rule fastani:
    input:
        #mags = expand("data/combined_bins/{MAG}.fasta",MAG = glob_wildcards('data/combined_bins/{MAG}.fasta').MAG),
        ref_genomes = expand("data/reference/reference_genomes/{genome}.fasta",genome = glob_wildcards('data/reference/reference_genomes/{genome}.fasta').genome)
    output:
        "data/fastANI/ani.txt"
    log:
        "logs/fastANI.log"
    conda: "code/fastani.yaml"
    resources: cpus=48, mem_mb=16000, time_min=120
    shell:
        """
        ls data/reference/reference_genomes/*.fasta > data/fastANI/genome_list.txt
        fastANI --ql data/fastANI/genome_list.txt --rl data/fastANI/genome_list.txt -t {resources.cpus} -o {output}
        """


rule phyloflash:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
    output: "data/phyloFlash/{sample}/{sample}.phyloFlash.html"
    params:
        database = "/geomicro/data2/kiledal/references/phyloflash/138.1",
        percent_id = 85,
        ref_clusterid = 100,
        taxlvl = 6
    conda: "code/phyloFlash.yaml"
    log: "logs/phyloFlash/{sample}_phyloFlash.log"
    benchmark: "benchmarks/phyloFlash/{sample}_phloFlash.txt"
    resources: cpus = 16, mem_mb = 20000, time_min = 1440
    shell:
        """
        mkdir -p data/phyloFlash/{wildcards.sample}
        cd data/phyloFlash/{wildcards.sample}
        
        phyloFlash.pl \
            -dbhome {params.database} \
            -lib {wildcards.sample} \
            -CPUs {resources.cpus} \
            -id {params.percent_id} \
            -clusterid {params.ref_clusterid} \
            -poscov \
            -almosteverything \
            -taxlevel {params.taxlvl}\
            -read1 {input.f_reads} \
            -read2 {input.r_reads} > ../../../{log} 2>&1

        tar -xzf *.phyloFlash.tar.gz
        """


checkpoint get_Microcystis_16S:
    input: 
        phyloflash_out = expand("data/phyloFlash/{sample}/{sample}.phyloFlash.html", sample = samples),
        script = "code/get_Phyloflash_Microcystis_sequences.R"
    output: 
        #expand("data/phyloFlash/{sample}/{sample}_microcystis_16S.fasta",sample = samples),
        #all_samples = 
        "data/phyloFlash/microcystis_seqs.fasta"
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        Rscript {input.script}
        """


def microcystis_16S_mapping(wildcards):
    checkpoint_output = checkpoints.get_Microcystis_16S.get(**wildcards).output[0]    
    file_names = expand("data/phyloFlash/read_mapping/{sample2}.sam", 
                        sample2= glob_wildcards("data/phyloFlash/{sample}/{sample2}_microcystis_16S.fasta").sample2)
    return file_names




rule map_reads_to_16S:
    input:
        "data/phyloFlash/microcystis_seqs.fasta",
        f_reads = "data/qcd_reads/{sample2}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample2}_R2_deduped_trimmed_screened.fq",
        ref = "data/phyloFlash/{sample2}/{sample2}_microcystis_16S.fasta"
    output:
        sam = "data/phyloFlash/read_mapping/{sample2}.sam"
    params: "idfilter=0.95 minidentity=0.9 nodisk=t po=f noheader=t saa=f nmtag=t idtag=t secondary=f"
    conda: "code/bbmap.yaml"
    log: "logs/phyloFlash/{sample2}_bbmap.log"
    resources: cpus=1
    shell:
        """
         #outm={output.sam}

         bbmap.sh ref={input.ref} in={input.f_reads} in2={input.r_reads} threads={resources.cpus} {params} outm=stdout.sam | reformat.sh in=stdin.sam out={output.sam} minidfilter=0.95 > {log} 2>&1
        """

rule run_phyloflash:
    input: 
        expand("data/phyloFlash/{sample}/{sample}.phyloFlash.html", sample = samples),
        #expand("data/phyloFlash/read_mapping/{sample}.sam",sample = samples)
        "data/phyloFlash/microcystis_seqs.fasta",
        microcystis_16S_mapping

# For variant calling
rule map_reads_to_marker_genes:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        ref = "data/reference/marker_genes/single_ref/{gene}.fa"
    output:
        int_sam = temp("data/marker_genes/{sample}/{gene}_mapped_temp.sam"),
        sam = "data/marker_genes/{sample}/{gene}_mapped.sam",
        bam = "data/marker_genes/{sample}/{gene}_mapped.bam",
        unsorted_bam = temp("data/marker_genes/{sample}/{gene}_mapped_unsorted.bam")
    params: "idfilter=0.95 minidentity=0.9 nodisk=t po=f noheader=f saa=f nmtag=t idtag=t secondary=f outputunmapped=f"
    conda: "code/bbmap.yaml"
    log: "logs/read_mapping/{sample}.{gene}_bbmap.log"
    resources: cpus=1
    shell:
        """
        bbmap.sh ref={input.ref} in={input.f_reads} in2={input.r_reads} threads={resources.cpus} {params} outm={output.int_sam} > {log} 2>&1

        # ID filter in bbmap didn't seem to be working, so filtering separately
        reformat.sh minidfilter=0.95 in={output.int_sam} out={output.sam} >> {log} 2>&1

        # If no headers: -bT {input.ref}
        samtools view -bS {output.sam} > {output.unsorted_bam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        samtools index {output.bam}
        """

# For spades assembly, only map paired reads
rule map_reads_to_marker_genes_paired_only:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        ref = "data/reference/marker_genes/multiple_refs/{gene}.fa"
    output:
        f_mapped = temp("data/marker_genes/{sample}/{gene}_f_mapped.fastq"),
        r_mapped = temp("data/marker_genes/{sample}/{gene}_r_mapped.fastq")
        #int_sam = temp("data/marker_genes/{sample}/{gene}_mapped_temp_paired.sam")
    params: "idfilter=0.95 minidentity=0.93 nodisk=t po=t noheader=t saa=f nmtag=t idtag=t secondary=f outputunmapped=f"
    conda: "code/bbmap.yaml"
    log: "logs/read_mapping/{sample}.{gene}_bbmap_PairedOnly.log"
    resources: cpus=1
    shell:
        """
        bbmap.sh ref={input.ref} in={input.f_reads} in2={input.r_reads} threads={resources.cpus} {params} out={output.f_mapped} out2={output.r_mapped}  > {log} 2>&1
        """
        #outm={output.int_sam}
        #reformat.sh minidfilter=0.95 in={output.int_sam} out={output.f_mapped} out2={output.r_mapped} >> {log} 2>&1

rule assemble_marker_genes:
    input: 
        f_mapped = rules.map_reads_to_marker_genes_paired_only.output.f_mapped,
        r_mapped = rules.map_reads_to_marker_genes_paired_only.output.r_mapped
    output:
        gene_assembly = "data/marker_genes/{sample}/{gene}_assembly.fasta"
    params:
        assembly_dir = "data/marker_genes/{sample}/{gene}_assembly",
        temp_contigs = "data/marker_genes/{sample}/{gene}_assembly/contigs.fasta"
    conda: "code/main_env.yaml"
    log: "logs/read_mapping/{sample}.{gene}_spades.log"
    resources: cpus=4
    shell:
        """
        spades.py -t {resources.cpus} -k 99,111,127 -1 {input.f_mapped} -2 {input.r_mapped} -o {params.assembly_dir} > {log} 2>&1

        cp {params.temp_contigs} {output.gene_assembly}; rm -r {params.assembly_dir}
        """

rule freebayes:
    input: 
        bam = rules.map_reads_to_marker_genes.output.bam,
        ref = rules.map_reads_to_marker_genes.input.ref
    output: "data/marker_genes/{sample}/{gene}.vcf"
    conda: "code/freebayes.yaml"
    log: "logs/read_mapping/{sample}.{gene}_freebayes.log"
    resources: cpus=1
    shell:
        """
        freebayes -f {input.ref} {input.bam} > {output}
        """

rule vcf_consensus:
    input: 
        ref = rules.map_reads_to_marker_genes.input.ref,
        vcf = rules.freebayes.output
    output: 
        compressed_vcf = "data/marker_genes/{sample}/{gene}.vcf.gz",
        consensus_fasta = "data/marker_genes/{sample}/{gene}_{sample}_consensus.fasta"
    conda: "code/vcftools.yaml"
    log: "logs/read_mapping/{sample}.{gene}_vcf_consensus.log"
    resources: cpus=1
    shell:
        """
        #bgzip -c {input.vcf} {output.compressed_vcf}
        
        bcftools view -Oz -o {output.compressed_vcf} {input.vcf}
        #htsfile {output.compressed_vcf}
        bcftools index {output.compressed_vcf}


        #tabix -p vcf {input.vcf}
        cat {input.ref} | vcf-consensus {output.compressed_vcf} > {output.consensus_fasta}

        sed -i "1s/.*/>{wildcards.gene}_{wildcards.sample}_consensus/" {output.consensus_fasta}
        """

rule tidy_fastas:
    input: expand("data/marker_genes/{sample}/{gene}_{sample}_consensus.fasta",sample = samples, gene = glob_wildcards("data/reference/marker_genes/single_ref/{gene}.fa").gene)
    output: expand("data/marker_genes/identified_variants/{gene_name}_called_variants.fasta", gene_name = glob_wildcards("data/reference/marker_genes/single_ref/{gene}.fa").gene)
    container: "docker://eandersk/r_microbiome"
    shell:
        """
        Rscript code/tidy_marker_genes.R
        """

rule remap_marker_genes_to_consensus:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        ref = rules.vcf_consensus.output.consensus_fasta
    output:
        int_sam = temp("data/marker_genes/{sample}/{gene}_mapped_to_consensus_temp.sam"),
        sam = "data/marker_genes/{sample}/{gene}_mapped_to_consensus.sam",
        bam = "data/marker_genes/{sample}/{gene}_mapped_to_consensus.bam",
        unsorted_bam = temp("data/marker_genes/{sample}/{gene}_mapped_to_consensus_unsorted.bam")
    params: "idfilter=0.95 minidentity=0.9 nodisk=t po=f noheader=f saa=f nmtag=t idtag=t secondary=f outputunmapped=f"
    conda: "code/bbmap.yaml"
    log: "logs/read_mapping/{sample}.{gene}_bbmap.log"
    resources: cpus=1
    shell:
        """
        bbmap.sh ref={input.ref} in={input.f_reads} in2={input.r_reads} threads={resources.cpus} {params} outm={output.int_sam} > {log} 2>&1

        # ID filter in bbmap didn't seem to be working, so filtering separately
        reformat.sh minidfilter=0.95 in={output.int_sam} out={output.sam} >> {log} 2>&1

        # If no headers: -bT {input.ref}
        samtools view -bS {output.sam} > {output.unsorted_bam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        samtools index {output.bam}
        """

rule remap_marker_genes_to_assembly:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        ref = rules.assemble_marker_genes.output.gene_assembly
    output:
        int_sam = temp("data/marker_genes/{sample}/{gene}_mapped_to_assembly_temp.sam"),
        sam = "data/marker_genes/{sample}/{gene}_mapped_to_assembly.sam",
        bam = "data/marker_genes/{sample}/{gene}_mapped_to_assembly.bam",
        unsorted_bam = temp("data/marker_genes/{sample}/{gene}_mapped_to_assembly_unsorted.bam")
    params: "idfilter=0.95 minidentity=0.9 nodisk=t po=f noheader=f saa=f nmtag=t idtag=t secondary=f outputunmapped=f"
    conda: "code/bbmap.yaml"
    log: "logs/read_mapping/{sample}.{gene}_bbmap.log"
    resources: cpus=1
    shell:
        """
        bbmap.sh ref={input.ref} in={input.f_reads} in2={input.r_reads} threads={resources.cpus} {params} outm={output.int_sam} > {log} 2>&1

        # ID filter in bbmap didn't seem to be working, so filtering separately
        reformat.sh minidfilter=0.95 in={output.int_sam} out={output.sam} >> {log} 2>&1

        # If no headers: -bT {input.ref}
        samtools view -bS {output.sam} > {output.unsorted_bam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        samtools index {output.bam}
        """

rule align_marker_variants:
    input: "data/marker_genes/identified_variants/{gene}_called_variants.fasta"
    output: "data/marker_genes/identified_variants/{gene}_called_variants.afa"
    log: "logs/align_variants_muscle/{gene}_muscle.log"
    conda: "code/muscle.yaml"
    resources: cpus=4
    shell:
        """
        muscle -align {input} -threads {resources.cpus} -output {output} > log
        """

rule calc_trees_from_variants:
    input: "data/marker_genes/identified_variants/{gene}_called_variants.afa"
    output: "data/marker_genes/identified_variants/{gene}_called_variants.tre"
    log: "logs/fasttree_variants/{gene}_fasttree.log"
    conda: "code/fasttree.yaml"
    shell:
        """
        FastTree -gtr -nt {input} > {output}
        """


#["cpcB","ftsZ","glnA","gltX","ntcA","pgi","recA","tpiA","16S", "itsC"]

rule do_read_mapping:
    input: 
        expand("data/marker_genes/{sample}/{gene}_f_mapped.fastq",sample = samples, gene = ["cpcB","ftsZ","glnA","gltX","ntcA","pgi","recA","tpiA","16S", "itsC"]),
        expand("data/marker_genes/{sample}/{gene}_assembly.fasta", sample = samples, gene = ["cpcB","ftsZ","glnA","gltX","ntcA","pgi","recA","tpiA","16S", "itsC"]),
        expand("data/marker_genes/{sample}/{gene}_{sample}_consensus.fasta", sample = samples, gene = glob_wildcards("data/reference/marker_genes/single_ref/{gene}.fa").gene),
        expand("data/marker_genes/{sample}/{gene}_mapped_to_assembly.sam", sample = samples, gene = ["cpcB","ftsZ","glnA","gltX","ntcA","pgi","recA","tpiA","16S", "itsC"]),
        expand("data/marker_genes/{sample}/{gene}_mapped_to_consensus.sam", sample = samples, gene = glob_wildcards("data/reference/marker_genes/single_ref/{gene}.fa").gene)

rule build_gene_trees:
    input:
        expand("data/marker_genes/identified_variants/{gene}_called_variants.tre",gene = glob_wildcards("data/reference/marker_genes/single_ref/{gene}.fa").gene)


rule gtotree:
    input:
        "data/reference/reference_genomes"
    output: 
        genome_list = "data/reference/reference_genomes/list.txt",
        gtotree_dir = directory("data/gtotree"),
        tree = "data/gtotree/gtotree.tre"
    conda: "code/gtotree.yaml"
    log: "logs/gtotree.log"
    benchmark: "benchmarks/gtotree.txt"
    resources: cpus=24
    shell:
        """
        ls {input}/*.fasta > {output.genome_list}

        GToTree -f {output.genome_list} \
            -H Cyanobacteria \
            -j {resources.cpus} \
            -F -o {output.gtotree_dir} > {log} 2>&1
        """


rule get_phylomark:
    output: 
        dir = directory("code/Phylomark"),
        phylomark = "code/Phylomark/phylomark.py"
    conda: "code/phylomark.yaml"
    shell:
        """
        rm -rf {output}
        cd code
        git clone https://github.com/jasonsahl/Phylomark.git
        cd Phylomark
        chmod +x phylomark.py
        python setup.py install --user
        """

rule phylomark_400:
    input:
        ref = "data/reference/reference_genomes/GCF_002095975.1_ASM209597v1_genomic.fasta",
        genomes = "data/reference/reference_genomes",
        tree = rules.gtotree.output.tree,
        phylomark = rules.get_phylomark.output.phylomark
    output:
        dir = directory("data/phylomark/400bp"),
        results =  "data/phylomark/400bp/results.txt",
        trees = "data/phylomark/400bp/all_trees.txt"
    conda: "code/phylomark.yaml"
    log: "logs/phylomark/400bp.log"
    benchmark: "benchmarks/phylomark/400bp.txt"
    params: frag_len = 400
    resources: cpus=80
    shell:
        """
        export PYTHONPATH=$PWD/code/Phyomark
        
        PROJ_DIR=$PWD
        PM=$PWD/code/Phylomark/phylomark.py

        cd $PWD/{output.dir}
        $PM \
            -r $PROJ_DIR/{input.ref} \
            -d $PROJ_DIR/{input.genomes} \
            -t $PROJ_DIR/{input.tree} \
            -l {params.frag_len} \
            -p {resources.cpus}
        """

rule phylomark_marker_genes:
    input:
        ref = "data/reference/reference_genomes/GCF_002095975.1_ASM209597v1_genomic.fasta",
        genomes = "data/reference/reference_genomes",
        tree = rules.gtotree.output.tree,
        phylomark = rules.get_phylomark.output.phylomark
    output:
        dir = directory("data/phylomark/marker_genes"),
        results =  "data/phylomark/marker_genes/results.txt",
        trees = "data/phylomark/marker_genes/all_trees.txt"
    conda: "code/phylomark.yaml"
    log: "logs/phylomark/marker_genes.log"
    benchmark: "benchmarks/phylomark/marker_genes.txt"
    params: frag_len = 400
    resources: cpus=80
    shell:
        """
        export PYTHONPATH=$PWD/code/Phyomark
        
        PROJ_DIR=$PWD
        PM=$PWD/code/Phylomark/phylomark.py

        cd $PWD/{output.dir}
        $PM \
            -g $PROJ_DIR/{input.ref} \
            -d $PROJ_DIR/{input.genomes} \
            -t $PROJ_DIR/{input.tree} \
            -p {resources.cpus}
        """


# # Unifrac distance calculations
# """
# biom convert --to-hdf5 -i data/bin_abund.biom -o data/bin_abund.hdf5.biom
# ssu -i bin_abund.hdf5.biom -t bin_ani_tree.nwk -m unweighted -o unweighted_unifrac.dst
# ssu -i bin_abund.hdf5.biom -t bin_ani_tree.nwk -m weighted_normalized -o weighted_unifrac.dst
# """


########################################
# BLAST searches for genes of interest #
########################################

rule run_BLAST:
    input:
        expand("data/BLAST/{blast_name}_nuc/{MAG}.blastn", MAG = glob_wildcards('data/combined_bins/{MAG}.fasta').MAG, blast_name = ["mlr","catalase", "katG"]),
        expand("data/BLAST/{blast_name}_prot/{MAG}.blastp", MAG = glob_wildcards('data/combined_bins/{MAG}.fasta').MAG, blast_name = ["mlr","t2prx"])


rule make_nuc_blastdbs:
    input:
        db = "data/reference/BLAST/nuc/{blast_name}.fasta",
    output:
        db_index = "data/reference/BLAST/nuc/{blast_name}.fasta.nin",
    output: directory("data/gtdb")
    log:
        "logs/BLAST/{blast_name}_nuc_makeblastdb.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        makeblastdb -in {input.db} -dbtype nucl -logfile logs/BLAST/{wildcards.blast_name}_nuc.dblog
        """

rule make_prot_blastdbs:
    input:
        db = "data/reference/BLAST/prot/{blast_name}.fasta"
    output:
        db_index = "data/reference/BLAST/prot/{blast_name}.fasta.pin"
    output: directory("data/gtdb")
    log:
        "logs/BLAST/{blast_name}_prot_makeblastdb.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        makeblastdb -in {input.db} -dbtype prot -logfile logs/BLAST/{wildcards.blast_name}_prot.dblog
        """

rule blast_nuc:
    input:
        blast_db = rules.make_nuc_blastdbs.input.db,
        blast_db_index = rules.make_nuc_blastdbs.output.db_index,
        mag = "data/combined_bins/{MAG}.fasta"
    output:
        blast_res = "data/BLAST/{blast_name}_nuc/{MAG}.blastn"
        #post_blast_res = "data/BLAST/{blast_name}_nuc/{MAG}.pb"
    log:
        "logs/BLAST/{blast_name}_nuc_{MAG}.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        blastn -query {input.mag} \
            -db {input.blast_db} \
            -out {output.blast_res} \
            -outfmt '6 std qcovs stitle qseq sseq' \
            -num_threads {resources.cpus}

        # removed additional output columns: -outfmt '6 std qcovs stitle' \
        # also removed database size standardization: -dbsize 1000000 \
        """
        #perl code/postBlast.pl -b {output.blast_res} -p 95 -s 50 -e 1e-5 -o {output.post_blast_res}

rule blast_prot:
    input:
        blast_db = rules.make_prot_blastdbs.input.db,
        blast_db_index = rules.make_prot_blastdbs.output.db_index,
        prodigal_mag = rules.prodigal_mags.params.faa
    output:
        blast_res = "data/BLAST/{blast_name}_prot/{MAG}.blastp"
        #post_blast_res = "data/BLAST/{blast_name}_prot/{MAG}.pb"
    log:
        "logs/BLAST/{blast_name}_prot_{MAG}.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        blastp -query {input.prodigal_mag} \
            -db {input.blast_db} \
            -out {output.blast_res} \
            -outfmt '6 std qcovs stitle qseq sseq' \
            -num_threads {resources.cpus}

        # removed additional output columns: -outfmt '6 std qcovs stitle' \
        # also removed database size standardization: -dbsize 1000000 \

        """
        #perl code/postBlast.pl -b {output.blast_res} -p 95 -s 50 -e 1e-5 -o {output.post_blast_res}


rule strainge_kmerize_genomes:
    input:
        "/geomicro/data2/kiledal/references/microcystis_genomes_from_jacob/{genome}.fa"
    output:
        "data/reference/strainge/reference_genomes/{genome}.hdf5"
    log: "logs/strainge_kmerize/{genome}.log"
    benchmark: "benchmarks/strainge_kmerize/{genome}.txt"
    conda: "code/strainge.yaml"
    shell:
        """
        straingst kmerize -o {output} {input}
        """

rule kmerize_genomes:
    input: expand("data/reference/strainge/reference_genomes/{genome}.hdf5", genome = glob_wildcards("/geomicro/data2/kiledal/references/microcystis_genomes_from_jacob/{genome}.fa").genome)


rule strainge_compare_refs:
    input:
        expand("data/reference/strainge/reference_genomes/{genome}.hdf5", genome = glob_wildcards("/geomicro/data2/kiledal/references/microcystis_genomes_from_jacob/{genome}.fa").genome)
    output:
        "data/strainge/ref_similarities.tsv"
    log: "logs/strainge_ref_compare.log"
    benchmark: "benchmarks/strainge_ref_compare.txt"
    conda: "code/strainge.yaml"
    shell:
        """
        straingst kmersim --all-vs-all -t 4 -S jaccard -S subset {input} > {output}
        """
    
rule strainge_cluster_refs:
    input:
        kmerized_genomes = expand("data/reference/strainge/reference_genomes/{genome}.hdf5", genome = glob_wildcards("/geomicro/data2/kiledal/references/microcystis_genomes_from_jacob/{genome}.fa").genome),
        ref_similarities = "data/strainge/ref_similarities.tsv"
    output:
        "data/strainge/references_to_keep.txt"
    log: "logs/strainge_ref_compare.log"
    benchmark: "benchmarks/strainge_ref_compare.txt"
    conda: "code/strainge.yaml"
    shell:
        """
        straingst cluster -i {input.ref_similarities} -d -C 0.99 -c 0.90 \
            --clusters-out clusters.tsv \
            {input.kmerized_genomes} > {output}
        """

rule strainge_make_pan_db:
    input:
        refs_to_use = "data/strainge/references_to_keep.txt"
    output:
        "data/strainge/pan-genome-db.hdf5"
    log: "logs/strainge_ref_compare.log"
    benchmark: "benchmarks/strainge_ref_compare.txt"
    conda: "code/strainge.yaml"
    shell:
        """
        straingst createdb -f {input.refs_to_use} -o {output}
        """

rule strainge_kmerize_reads:
    input:
        f_seq = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_seq = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq"
    output:
        "data/strainge/kmerized_samples/{sample}.hdf5"
    log: "logs/strainge_kmerize_reads/{sample}.log"
    benchmark: "benchmarks/strainge_kmerize_reads/{sample}.txt"
    conda: "code/strainge.yaml"
    shell:
        """
        straingst kmerize -k 23 -o {output} \
            {input.f_seq} {input.r_seq}
        """

rule run_read_kmerization:
    input: expand("data/strainge/kmerized_samples/{sample}.hdf5", sample = samples)

rule strainge_run:
    input:
        db = "data/strainge/pan-genome-db.hdf5",
        sample_kmers = "data/strainge/kmerized_samples/{sample}.hdf5"
    params:
        out_prefix = "data/strainge/results/{sample}"
    output:
        strains = "data/strainge/results/{sample}.strains.tsv",
        stats = "data/strainge/results/{sample}.stats.tsv"
    resources: mem_mb = 50000, cpus = 1
    log: "logs/strainge/{sample}.log"
    benchmark: "benchmarks/strainge/{sample}.txt"
    conda: "code/strainge.yaml"
    shell:
        """
        straingst run -O -o {params.out_prefix} {input.db} {input.sample_kmers}
        """



rule strainge_kmerize_2014_metaG_reads:
    input:
        f_seq = "/geomicro/data2/kiledal/GLAMR/import/metagenomes/{sample}_1.fastq.gz",
        r_seq = "/geomicro/data2/kiledal/GLAMR/import/metagenomes/{sample}_2.fastq.gz"
    output: "data/strainge/kmerized_samples/{sample}.hdf5"
    log: "logs/strainge_kmerize_reads/{sample}.log"
    benchmark: "benchmarks/strainge_kmerize_reads/{sample}.txt",
    resources: mem_mb = 75000, cpus = 1, time_min = 1000
    conda: "code/strainge.yaml"
    shell:
        """
        export HDF5_USE_FILE_LOCKING='FALSE' # prevents lock error on NFS when writing output

        straingst kmerize -k 23 -o {output} \
            {input.f_seq} {input.r_seq} 1>{log} 2>&1
        """

rule run_strainge:
    input:
        expand("data/strainge/results/{sample}.strains.tsv", sample = samples),
        expand("data/strainge/results/{metag}.strains.tsv", metag = glob_wildcards("/geomicro/data2/kiledal/GLAMR/import/metagenomes/{sample}_1.fastq.gz").sample)



###### Marker gene analysis for Laura ####


rule make_nuc_blastdb_markers:
    input: 
        genome = "data/reference/laura_ref_genomes/{genome}.fna"
    output:
        #concat_genomes = "data/reference/laura_ref_genomes/combined_db.fasta",
        db_index = "data/reference/laura_ref_genomes/{genome}.fna.nin"
    log:
        "logs/BLAST/laura_genomes_nuc_makeblastdb_{genome}.log"
    resources: cpus=1, mem_mb=5000, time_min=120
    shell:
        """
        makeblastdb -in {input.genome} -dbtype nucl -logfile logs/BLAST/laura_genomes_db_nuc_{wildcards.genome}.dblog
        """
         #cat data/reference/laura_ref_genomes/*.fna > {output.concat_genomes}

rule blast_nuc_laura:
    input:
        blast_db = rules.make_nuc_blastdb_markers.input.genome,
        blast_db_index = rules.make_nuc_blastdb_markers.output.db_index,
        gene = "data/reference/marker_genes/multiple_refs/{gene}.fa"
    output:
        blast_res = "data/BLAST/laura_markers/{gene}__{genome}.blastn"
        #post_blast_res = "data/BLAST/{blast_name}_nuc/{MAG}.pb"
    log:
        "logs/BLAST/laura_markers/{gene}__{genome}.log"
    resources: cpus=8, mem_mb=5000, time_min=120
    shell:
        """
        blastn -query {input.gene} \
            -db {input.blast_db} \
            -out {output.blast_res} \
            -outfmt '6 std qcovs stitle qseq sseq' \
            -num_threads {resources.cpus}

        # removed additional output columns: -outfmt '6 std qcovs stitle' \
        # also removed database size standardization: -dbsize 1000000 \
        """

rule run_laura_BLASTS:
    input: expand("data/BLAST/laura_markers/{gene}__{genome}.blastn", gene = glob_wildcards("data/reference/marker_genes/multiple_refs/{gene}.fa").gene, genome = glob_wildcards("data/reference/laura_ref_genomes/{genome}.fna").genome)




rule mlr_hmm:
    input:
        mag_genes = "data/prodigal_mag/{MAG}.faa",
        hmms = "data/reference/hmms/{gene}.hmm"
    output:
        hits = "data/hmms/{gene}/{MAG}.txt",
        hits_preprocessed = "data/hmms/{gene}/{MAG}__formatted.txt"
    log:
        "logs/hmmsearch/{gene}/{MAG}.log"
    resources: cpus=8, mem_mb=5000, time_min=120
    shell:
        """
        hmmsearch --cpu {resources.cpus} --tblout {output.hits} {input.hmms} {input.mag_genes} > {log} 2>&1

        tail -n+4 {output.hits} | head -n -10 | sed 's/ * / /g' | cut -f 1-10 -d " " > {output.hits_preprocessed}
        """


rule run_hmmsearch:
    input: 
        expand("data/hmms/{gene}/{MAG}.txt", gene = "mlr", MAG = glob_wildcards("data/combined_bins/{MAG}.fasta").MAG),
        expand("data/hmms/{gene}/{MAG}__formatted.txt", gene = "mlr", MAG = glob_wildcards("data/combined_bins/{MAG}.fasta").MAG)


rule kofam_scan:
    input:
        genes = "data/prodigal_mag/{MAG}.renamed.faa",
        profile = "../../GLAMR/data/reference/kegg/kofamscan/profiles",
        ko_list = "../../GLAMR/data/reference/kegg/kofamscan/ko_list"
    output:
        ko_annot = "data/kofam_scan/{MAG}__kofam_results.txt"
    conda: "code/kofamscan.yaml"
    #shadow: "shallow"
    benchmark: "benchmarks/kofamscan/{MAG}.txt"
    log: "logs/kofamscan/{MAG}.log"
    resources: cpus=24, time_min = 20000, mem_mb = lambda wildcards, attempt: attempt * 20000
    shell:
        """
        exec_annotation \
            -o {output.ko_annot} \
            --cpu={resources.cpus}  \
            -f mapper \
            --profile {input.profile} \
            --tmp-dir=/tmp/{wildcards.MAG}_kofamscan \
            --ko-list {input.ko_list} {input.genes}
        """
        

rule run_kofam_scan:
    input:
        expand("data/kofam_scan/{MAG}__kofam_results.txt", MAG = glob_wildcards("data/combined_bins/{MAG}.fasta").MAG)

rule KEGGdecoder:
    input: "data/kofam_scan/{MAG}__kofam_results.txt"
    output: "data/kegg_decoder/{MAG}/kegg_decoder_list.tsv"
    conda: "code/kegg_decoder.yml"
    shell:
        """
        KEGG-decoder --input {input} --output {output} --vizoption static        
        """

rule run_KEGGdecoder:
    input:
        expand("data/kegg_decoder/{MAG}/kegg_decoder_list.tsv", MAG = glob_wildcards("data/combined_bins/{MAG}.fasta").MAG)



rule rename_fasta:
    input: 
        script = "code/rename_megahit_contigs.R",
        #assembly_dir = "data/omics/metagenomes/{sample}/assembly/megahit",
        contigs = "data/prodigal_mag/{MAG}.faa",
        #assembly_done = "data/omics/metagenomes/{sample}/assembly/megahit/.done"
    output:
        contigs = "data/prodigal_mag/{MAG}.renamed.faa",
        #done = touch("data/omics/metagenomes/{sample}/assembly/megahit/.contigs_renamed")
    params:
        work_dir = "/geomicro/data2/kiledal/projects/microcystis_genomes_2022"
    container: "docker://eandersk/r_microbiome"
    benchmark: "benchmarks/rename_megahit_contigs/{MAG}.txt"
    resources: cpus = 1, time_min=2000, mem_mb = 75000
    shell:
        """
        cd {params.work_dir}
        
        pwd
        
        ./{input.script} {input.contigs}
        """



# For gene P/A determination
rule map_reads_to_genes:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        ref = "data/reference/read_mapping_refs/{genePA}.fa"
    output:
        int_sam = temp("data/gene_PA_read_mapping/{sample}/{genePA}_mapped_temp.sam"),
        sam = "data/gene_PA_read_mapping/{sample}/{genePA}_mapped.sam",
        bam = "data/gene_PA_read_mapping/{sample}/{genePA}_mapped.bam",
        unsorted_bam = temp("data/gene_PA_read_mapping/{sample}/{genePA}_mapped_unsorted.bam")
    params: "idfilter=0.95 minidentity=0.9 nodisk=t po=f noheader=f saa=f nmtag=t idtag=t secondary=f outputunmapped=f"
    conda: "code/bbmap.yaml"
    log: "logs/read_mapping/{sample}.{genePA}_bbmap.log"
    resources: cpus=1
    shell:
        """
        bbmap.sh ref={input.ref} in={input.f_reads} in2={input.r_reads} threads={resources.cpus} {params} outm={output.int_sam} > {log} 2>&1

        # ID filter in bbmap didn't seem to be working, so filtering separately
        reformat.sh minidfilter=0.95 in={output.int_sam} out={output.sam} >> {log} 2>&1

        # If no headers: -bT {input.ref}
        samtools view -bS {output.sam} > {output.unsorted_bam}
        samtools sort -o {output.bam} {output.unsorted_bam}
        samtools index {output.bam}
        """

rule run_map_reads_to_genes:
    input: expand("data/gene_PA_read_mapping/{sample}/{genePA}_mapped.bam", sample = samples, genePA = ["mcyGenes"])

rule map_reads_to_genes2:
    input:
        f_reads = "data/qcd_reads/{sample}_R1_deduped_trimmed_screened.fq",
        r_reads = "data/qcd_reads/{sample}_R2_deduped_trimmed_screened.fq",
        ref = "data/reference/read_mapping_refs/{genePA}.fa"
    output:
        temp_bam = temp("data/gene_PA_read_mapping_BWA/{sample}/{genePA}_mapped_temp.bam"),
        sam = temp("data/gene_PA_read_mapping_BWA/{sample}/{genePA}_mapped.sam"),
        bam = "data/gene_PA_read_mapping_BWA/{sample}/{genePA}_mapped.bam",
        unsorted_bam = temp("data/gene_PA_read_mapping_BWA/{sample}/{genePA}_mapped_unsorted.bam")
    conda: "code/bwa.yaml"
    log: "logs/read_mapping/{sample}.{genePA}_BWA.log"
    resources: cpus=8
    shell:
        """
        bwa-mem2 index {input.ref}
        
        bwa-mem2 mem \
            -t {resources.cpus} \
            {input.ref} \
            {input.f_reads} {input.r_reads} > {output.sam}

        samtools view -bS {output.sam} > {output.temp_bam}
        
        filterBam \
            --in {output.temp_bam} \
            --out {output.unsorted_bam} \
            --minCover 80 \
            --minId 95
        
        samtools sort -o {output.bam} -@ {resources.cpus} {output.unsorted_bam}
        samtools index -@ {resources.cpus} {output.bam}
        """

rule run_map_reads_to_genes_BWA:
    input: expand("data/gene_PA_read_mapping_BWA/{sample}/{genePA}_mapped.bam", sample = samples, genePA = ["mergedRefs"])

