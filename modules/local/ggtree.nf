#!/usr/bin/env nextflow

process ggtree {
    publishDir "${params.out_dir}/ggtree/", mode: "copy"
    label "low"

    input:
    path (treefiles)

    output:
    path "*.png",emit:png

    script:
    """
    for treefile in ${treefiles}; do
        filename=\$(basename "\$treefile")
        if [[ "\$filename" != *no_msa* ]]; then
            plot_tree.R "\$treefile"
        else 
            echo "no virus sequences found for this sample, skipping ggtree analysis" > "\$filename.png"
        fi
    done
    """
}
