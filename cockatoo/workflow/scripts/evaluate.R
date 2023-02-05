##################
### evaluate.R ###
##################
# Author: Samuel Aroney

library(gt)
library(cowplot)
library(tidyverse)

matched_hits <- read_tsv(snakemake@input[["matched_hits"]])
novel_hits <- read_tsv(snakemake@input[["novel_hits"]])
coassemble_summary <- read_tsv(snakemake@params[["coassemble_summary"]])

################
### Plotting ###
################
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

#####################
### Summary stats ###
#####################
recovered_hits <- bind_rows(
    novel_hits,
    matched_hits %>% filter(!is.na(genome))
)

stat_groups <- tribble(
    ~statistic, ~within,
    "sequences", "targets",
    "taxonomy", "targets",
    "nontarget_sequences", "recovery",
    "novel_sequences", "recovery",
    "bins", "recovery"
)

summary_stats <- analysis %>%
    group_by(coassembly, status = c("match", "nonmatch")[as.integer(is.na(genome)) + 1]) %>%
    summarise(
        sequences = n(),
        bins = sum(!is.na(unique(genome))),
        taxonomy = sum(!is.na(unique(taxonomy))),
    ) %>%
    left_join(
        recovered_hits %>%
            group_by(coassembly, status = c("match", "nonmatch")[as.integer(is.na(found_in)) + 1]) %>%
            summarise(nontarget_sequences = n())
        ) %>%
    left_join(
        recovered_hits %>%
            group_by(coassembly, status = c("match", "nonmatch")[as.integer(!is.na(found_in) | !is.na(target)) + 1]) %>%
            summarise(novel_sequences = n())
        ) %>%
    pivot_longer(-c(coassembly, status), names_to = "statistic", values_to = "value") %>%
    {
        .["value"][.["statistic"] == "bins" & .["status"] == "nonmatch"] <- length(snakemake@config$recovered_bins) - .["value"][.["statistic"] == "bins" & .["status"] == "match"]
        .
    } %>%
    pivot_wider(names_from = status, values_from = value) %>%
    mutate(
        total = match + nonmatch,
        match_percent = round(match / total * 100, 2)
    ) %>%
    left_join(stat_groups)

summary_stats %>%
    select(coassembly, statistic, within, match, nonmatch, total, match_percent) %>%
    mutate(statistic = factor(statistic, levels = c("sequences", "taxonomy", "nontarget_sequences", "novel_sequences", "bins"))) %>%
    arrange(coassembly, statistic) %>%
    write_tsv(snakemake@output[["summary_stats"]])

####################
### Table output ###
####################
summary_table <- analysis %>%
    gt()

gtsave(summary_table, snakemake@output[["summary_table"]])

# Save R image for further processing
save.image(file = str_c(main_dir, "/evaluate.RData"))
