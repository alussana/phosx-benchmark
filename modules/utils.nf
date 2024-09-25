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


/*
split a single input file into chunks of n rows,
each named with a unique identifier

typically, each chunk in the output channel is mapped to its id with
.flatten().map{ file -> tuple( file.baseName, file ) }
*/
process split {

    input:
        path 'input/file'
        val n

    output:
        path "output/*"

    script:
    """
    mkdir -p output
    
    split -l ${n} -a 16 -x input/file output/
    """

}


/*
concatenate input files removing empty lines
*/
process concatenate {

    input:
        path 'input/*'
    
    output:
        path 'output/file'
    
    script:
    """
    mkdir -p output

    cat input/* | awk 'NF' > output/file
    """

}


/*
translate all the words in the second field of input/file.tsv 
specified in the first tab-separated column of input/dict.tsv
with the corresponding word found in the second column

don't discard untranslated rows
*/
process translate2col {

    input:
        path 'input/file.tsv'
        path 'input/dict.tsv'

    output:
        path 'translated_file.tsv'

    script:
    """
    translator.py \
        input/dict.tsv \
        input/file.tsv \
        1 \
        3 \
        2 \
        1 \
        > translated_file.tsv
    """

}