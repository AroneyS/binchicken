##################
### evaluate.R ###
##################
# Author: Samuel Aroney

library(cowplot)
library(tidyverse)

matched_hits <- read_tsv(snakemake@input[["matched_hits"]])
novel_hits <- read_tsv(snakemake@input[["novel_hits"]])

# Plotting
main_dir <- snakemake@output[["plots_dir"]]
dir.create(main_dir, recursive = TRUE)
taxonomy_groups <- c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

analysis <- matched_hits %>%
    filter(!is.na(target)) %>%
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
recovered_hits <- bind_rows(
    novel_hits,
    matched_hits %>% filter(!is.na(genome))
)

summary_stats <- analysis %>%
    group_by(coassembly, status = c("recovered", "missed")[as.integer(is.na(genome)) + 1]) %>%
    summarise(
        sequence_targets = n(),
        bins = sum(!is.na(unique(genome))),
        taxonomy = sum(!is.na(unique(taxonomy))),
    ) %>%
    left_join(
        recovered_hits %>%
            group_by(coassembly, status = c("recovered", "missed")[as.integer(is.na(found_in)) + 1]) %>%
            summarise(nontarget_recovered = n())
        ) %>%
    left_join(
        recovered_hits %>%
            group_by(coassembly, status = c("recovered", "missed")[as.integer(!is.na(found_in) | !is.na(target)) + 1]) %>%
            summarise(novel_recovered = n())
        ) %>%
    pivot_longer(-c(coassembly, status), names_to = "statistic", values_to = "value") %>%
    pivot_wider(names_from = status, values_from = value) %>%
    mutate(
        total = recovered + missed,
        recovered_percent = round(recovered / total * 100, 2)
    )

summary_stats %>%
    mutate(statistic = factor(statistic, levels = c("sequence_targets", "nontarget_recovered", "novel_recovered", "bins", "taxonomy"))) %>%
    arrange(coassembly, statistic) %>%
    write_tsv(snakemake@output[["summary_stats"]])

# Save R image for further processing
save.image(file = str_c(main_dir, "/evaluate.RData"))
