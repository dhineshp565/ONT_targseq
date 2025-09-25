#!/usr/bin/env nextflow

process htmltopdf {
    publishDir "${params.out_dir}/htmltopdf/", mode: "copy"
    label "medium"

    input:
    path html

    output:
    path "output.pdf"

    script:
    """
    weasyprint pdf.html output.pdf
    """
}
