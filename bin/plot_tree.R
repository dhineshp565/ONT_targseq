#!/usr/bin/env Rscript
# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)
treefile <- args[1]

# Load required libraries
library(ggtree)
library(ape)
library(dplyr)
library(ggplot2)

# Read the phylogenetic tree from the input file
tree <- read.tree(treefile)

# Create a data frame for tip labels and assign colors based on label content
tip_df <- data.frame(label = tree$tip.label) %>%
    mutate(color = ifelse(grepl("PDS", label), "red", "black"))

# Generate the tree plot with customizations
p <- ggtree(tree, layout='rectangular', size=0.2) %<+% tip_df +
    geom_tiplab(aes(color = color), size = 2.0) +                # Color tip labels
    scale_color_identity() +                                     # Use specified colors directly
    xlim(0,1.5) +                                                  # Set x-axis limits
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size =1.5) + # Label internal nodes
    geom_treescale(x = 7, y = 10, width = 0.5, fontsize = 1.5, linesize = 0.2) + # Add tree scale
    theme_tree2() +                                              # Apply tree theme
    labs(caption="nucleotide substitutions per site") +          # Add caption
    theme(plot.caption = element_text(hjust = 0.5))              # Center caption

# Construct output file name based on input file
out_file <- paste0("tree_", tools::file_path_sans_ext(basename(treefile)), ".png")

# Save the plot to a PNG file
ggsave(out_file, p, width = 15, height = 10, units = "cm", dpi = 300)
