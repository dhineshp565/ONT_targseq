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

max_branch_length <- max (node.depth.edgelength(tree))
x_max <- max_branch_length * 1.1

# Generate the tree plot with customizations
p <- ggtree(tree, layout='rectangular', size=0.2) %<+% tip_df +
    geom_tiplab(aes(color = color), size = 2.5,hjust=0) +                # Color tip labels
    scale_color_identity() +                                     # Use specified colors directly
    #xlim(0,9.0) +                                                  # Set x-axis limits
    geom_text2(aes(subset = !isTip, label = label), hjust = -0.2, size =1.9) + # Label internal nodes
    #geom_treescale(x = 10, y = 15, width = 0.5, fontsize = 1.5, linesize = 0.2) + # Add tree scale
    #theme_tree2() +                                              # Apply tree theme
    #labs(caption="nucleotide substitutions per site") +  
    geom_treescale(x = 0, y = -1, width = 0.5, fontsize = 3, linesize = 0.5) +   
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +   
    coord_cartesian (xlim=c(0,x_max))+   # Add caption
    theme(plot.caption = element_text(hjust = 0.5)) 

# Construct output file name based on input file
out_file <- paste0(tools::file_path_sans_ext(basename(treefile)), ".png")

# Save the plot to a PNG file
ggsave(out_file, p, width = 20, height = 16, units = "cm", dpi = 300)