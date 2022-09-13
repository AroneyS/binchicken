##################
### evaluate.R ###
##################
# Author: Samuel Aroney

library(cowplot)
library(tidyverse)

unbinned_hits <- read_tsv(snakemake@input[["unbinned_hits"]])

# Plotting
main_dir <- snakemake@output[["plots_dir"]]
dir.create(main_dir, recursive = TRUE)
taxonomy_groups <- c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

analysis <- unbinned_hits %>%
    separate(taxonomy, into = taxonomy_groups, sep = "; ", remove = FALSE)

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
}

plot_bars(analysis, str_c(main_dir, "/combined"))

plot_coassembly <- function(coassembly) {
    analysis %>%
        filter(coassembly == {{coassembly}}) %>%
        plot_bars(str_c(main_dir, coassembly, sep = "/"))
}

analysis %>%
    distinct(coassembly) %>%
    `$`(coassembly) %>%
    lapply(plot_coassembly)

# Summary stats
analysis %>%
    group_by(coassembly, status = c("recovered", "missed")[as.integer(is.na(genome)) + 1]) %>%
    summarise(
        sequences = n(),
        bins = sum(!is.na(unique(genome))),
        taxonomy = sum(!is.na(unique(taxonomy))),
    ) %>%
    pivot_longer(-c(coassembly, status), names_to = "statistic", values_to = "value") %>%
    pivot_wider(names_from = status, values_from = value) %>%
    mutate(
        total = recovered + missed,
        recovered_percent = round(recovered / total * 100, 2)
    ) %>%
    write_tsv(snakemake@output[["summary_stats"]])
