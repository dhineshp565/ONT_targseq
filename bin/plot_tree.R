#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
treefile <- args[1]

library(ggtree)
library(ape)
library(dplyr)
library(ggplot2)

tree <- read.tree(treefile)

tip_df <- data.frame(label = tree$tip.label) %>%
    mutate(color = ifelse(grepl("PDS", label), "red", "black"))

p <- ggtree(tree, layout='rectangular', size=0.2) %<+% tip_df +
    geom_tiplab(aes(color = color), size = 2.0) +
    scale_color_identity() +
    xlim(0,9) +
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size =1.5) +
    geom_treescale(x = 7, y = 10, width = 0.5, fontsize = 1.5, linesize = 0.2) +
    theme_tree2() +
    labs(caption="nucleotide substitutions per site") +
    theme(plot.caption = element_text(hjust = 0.5))

out_file <- paste0("tree_", tools::file_path_sans_ext(basename(treefile)), ".png")
ggsave(out_file, p, width = 15, height = 10, units = "cm", dpi = 300)
