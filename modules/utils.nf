#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
publish a file giving it an arbitrary path
*/
process publish {

    publishDir "${out_dir}",
                pattern: "${outputFile}",
                mode: 'copy'

    input:
        tuple val(outputFile),
              file('input/file')

    output:
        path "${outputFile}"

    script:
    """
    mkdir -p \$(dirname ${outputFile})
    
    cp input/file ${outputFile}
    """

}