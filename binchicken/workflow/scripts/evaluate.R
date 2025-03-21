##################
### evaluate.R ###
##################
# Author: Samuel Aroney

test <- !is.null(snakemake@config$test)
if (test) sink(file("/dev/null", "w"), type = "message")

Sys.setenv(OPENSSL_CONF = "/dev/null")
library(gt)
library(cowplot)
library(tidyverse)
if (!webshot::is_phantomjs_installed()) webshot::install_phantomjs()

# Snakemake inputs
matched_hits_path <- snakemake@input[["matched_hits"]]
novel_hits_path <- snakemake@input[["novel_hits"]]
cluster_summary_path <- snakemake@input[["cluster_summary"]]
summary_stats_path <- snakemake@input[["summary_stats"]]

coassemble_summary_path <- snakemake@params[["coassemble_summary"]]

main_dir <- snakemake@output[["plots_dir"]]
summary_table_path <- snakemake@output[["summary_table"]]

# Load inputs
matched_hits <- read_tsv(matched_hits_path)
novel_hits <- read_tsv(novel_hits_path)
coassemble_summary <- read_tsv(coassemble_summary_path)
summary_stats <- read_tsv(summary_stats_path)

if (is.null(cluster_summary_path)) {
    cluster_summary <- tibble(type = c("original", coassemble_summary$coassembly), clusters = NA_real_)
} else {
    cluster_summary <- read_csv(cluster_summary_path, col_names = c("type", "clusters"))
}

#################
### Functions ###
#################
plot_bars <- function(df, output_dir) {
    dir.create(output_dir, recursive = TRUE)

    bar_layers <- list(
        geom_col(),
        scale_fill_brewer(breaks = c("recovered", "missed"), palette = "Set2"),
        ylab("SingleM sequences"),
        coord_flip(),
        theme_cowplot(),
        theme(axis.title.y = element_blank())
    )

    df %>%
        group_by(Phylum) %>%
        summarise(
            recovered = sum(!is.na(genome)),
            missed = sum(is.na(genome))
        ) %>%
        pivot_longer(-Phylum, names_to = "outcome", values_to = "count") %>%
        ggplot(aes(reorder(Phylum, desc(count), FUN = sum), count, fill = outcome)) +
        bar_layers
    ggsave(str_c(output_dir, "/phylum_recovered.png"), dpi = 900, height = 8, width = 13)

    df %>%
        group_by(gene) %>%
        summarise(
            recovered = sum(!is.na(genome)),
            missed = sum(is.na(genome))
        ) %>%
        pivot_longer(-gene, names_to = "outcome", values_to = "count") %>%
        mutate(gene_number = map_int(gene, ~ str_split(., "\\.")[[1]][2] %>% as.integer())) %>%
        ggplot(aes(reorder(gene, gene_number, FUN = sum), count, fill = outcome)) +
        bar_layers
    ggsave(str_c(output_dir, "/gene_recovered.png"), dpi = 900, height = 8, width = 13)

    return(invisible(NULL))
}

plot_coassembly <- function(coassembly) {
    analysis %>%
        filter(coassembly == {{coassembly}}) %>%
        plot_bars(str_c(main_dir, coassembly, sep = "/"))

    return(invisible(NULL))
}

################
### Plotting ###
################
dir.create(main_dir, recursive = TRUE)
taxonomy_groups <- c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

analysis <- matched_hits %>%
    filter(!is.na(target)) %>%
    separate(taxonomy, into = taxonomy_groups, sep = "; ", remove = FALSE)

if (test) warning(analysis)

plot_bars(analysis, str_c(main_dir, "/combined"))

coassembly_plots <- analysis %>%
    distinct(coassembly) %>%
    pull(coassembly) %>%
    map(plot_coassembly)

####################
### Table output ###
####################
target_summary <- summary_stats %>%
    filter(statistic != "taxonomy") %>%
    pivot_wider(id_cols = coassembly, names_from = statistic, values_from = match)

target_totals <- summary_stats %>%
    filter(!statistic %in% c("nontarget_unbin_sequences", "novel_sequences", "taxonomy")) %>%
    pivot_wider(id_cols = coassembly, names_from = statistic, values_from = total) %>%
    rename(total_targets = sequences, total_bins = bins, total_recovered = nontarget_bin_sequences)

target_percentage <- summary_stats %>%
    filter(!statistic %in% c("bins", "nontarget_bin_sequences", "nontarget_unbin_sequences", "novel_sequences", "taxonomy")) %>%
    pivot_wider(id_cols = coassembly, names_from = statistic, values_from = match_percent) %>%
    rename(perc_targets = sequences)

original_clusters <- cluster_summary %>%
    filter(type == "original") %>%
    pull(clusters)

cluster_info <- cluster_summary %>%
    pivot_wider(names_from = type, values_from = clusters) %>%
    pivot_longer(-original, names_to = "coassembly", values_to = "clusters") %>%
    mutate(novel_clusters = clusters - original) %>%
    select(coassembly, novel_clusters)

if (!"unmapped_size" %in% colnames(coassemble_summary)) coassemble_summary$unmapped_size <- NA_real_

summary_table <- coassemble_summary %>%
    filter(coassembly %in% analysis$coassembly) %>%
    select(-c(samples, total_targets)) %>%
    left_join(target_summary) %>%
    left_join(target_totals) %>%
    left_join(target_percentage) %>%
    left_join(cluster_info) %>%
    mutate(
        sequences = pmap_chr(list(sequences, total_targets, perc_targets), ~ str_c(..1, "/", ..2, " (", ..3, "%)")),
        bins = map2_chr(bins, total_bins, ~ str_c(.x, " (", .y, " total)")),
        total_size = total_size / 10**9,
        unmapped_size = unmapped_size / 10**9,
        ) %>%
    select(coassembly, length, total_size, unmapped_size, bins, sequences, nontarget_bin_sequences, nontarget_unbin_sequences, novel_sequences, novel_clusters) %>%
    gt() %>%
    tab_spanner(
        label = "Gbp",
        columns = c(total_size, unmapped_size)
    ) %>%
    tab_spanner(
        label = "Recovered sequences",
        columns = c(sequences, nontarget_bin_sequences, nontarget_unbin_sequences, novel_sequences)
    ) %>%
    fmt_integer(c(length, total_size, unmapped_size, nontarget_bin_sequences, nontarget_unbin_sequences, novel_sequences, novel_clusters)) %>%
    cols_label(
        coassembly = "coassembly",
        length = "samples",
        total_size = "size",
        unmapped_size = "update",
        bins = "target bins",
        sequences = "targets",
        nontarget_unbin_sequences = "non-targets",
        nontarget_bin_sequences = "prior binned",
        novel_sequences = "novel",
        novel_clusters = "novel clusters"
    ) %>%
    tab_header(
        title = "Bin Chicken coassembly evaluation",
        subtitle = str_c("Original genomes had ", original_clusters, " clusters")
        ) %>%
    tab_style(
        style = cell_fill(),
        locations = cells_body(columns = c(length, unmapped_size, sequences, novel_sequences))
        )

if (is.null(cluster_summary_path)) summary_table <- summary_table %>% cols_hide(novel_clusters)

gtsave(summary_table, summary_table_path)

# Save R image for further processing
save.image(file = str_c(main_dir, "/evaluate.RData"))
