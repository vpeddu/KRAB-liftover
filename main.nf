OUTDIR = "${params.OUTDIR}"
BEDFILE = file("${params.BED}")
//CHAINFILE = file ("${params.CHAIN}")
//REF_FASTA = file("${params.REFERENCE}")
//BASE = "${params.base}"


METADATA = file("${params.METADATA}")


metadataCh = Channel
        .fromPath(METADATA)
        .splitCsv(header:true)
        .map{ row-> tuple(row.blank, (row.SampleName), file(row.ChainFile), file(row.ReferenceFasta)) }

//metadataCh.view()  


process liftOver {  
    input: 
        file BEDFILE 
        tuple(val(blank), val(BASE), file(CHAINFILE), file(REFERENCE)) from metadataCh
        //file CHAINFILE
        //val BASE
    container 'quay.io/biocontainers/ucsc-liftover:357--h446ed27_4'   
    publishDir "${OUTDIR}/liftOver/${BASE}.liftover", mode: 'copy'

    output:
        file "${BASE}.transferred.bed" into liftoverOutcH
        file "${BASE}.unmapped.bed"
        tuple( val(BASE), file(CHAINFILE), file(REFERENCE)) into sampleCh
    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    liftOver ${BEDFILE} ${CHAINFILE} ${BASE}.transferred.bed ${BASE}.unmapped.bed

    """
}

//sampleCh.view()


process getFasta {  
    input: 
        file BEDFILE from liftoverOutcH
        tuple(val(BASE), file(CHAINFILE), file(REFERENCE)) from sampleCh
    container 'quay.io/biocontainers/bedtools:2.26.0gx--he513fc3_4'   
    publishDir "${OUTDIR}/getFasta/${BASE}.getFasta", mode: 'copy'

    output:
        tuple(val(BASE), file("${BASE}.extracted.fasta")) into toSmashCh

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    bedtools getfasta -s -name -fi ${REFERENCE} -bed ${BEDFILE} > ${BASE}.extracted.fasta


    """
}

process smashFasta {  
    input: 
        tuple(val(BASE),file(EXONFASTA)) from toSmashCh
        //val BASE

    container 'quay.io/vpeddu/rgeneratesummary:latest'   
    publishDir "${OUTDIR}/smashFasta", mode: 'copy'

    output:
        file "${BASE}.smashed.fasta"

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    Rscript --vanilla ${baseDir}/bin/smash_fasta.r ${EXONFASTA}

    mv smashed.fasta ${BASE}.smashed.fasta

    """
}


def helpMessage() {
    log.info"""
    KRAB liftover workflow 
    Usage:
    An example command for running the pipeline is as follows:
    nextflow run vpeddu/KRAB-liftover \\
        --BASE <basename>\\
        --CHAIN <chainfile> \\
        --BED <bedfile> \\
        --REF_FASTA <reference_fasta> \\
        --OUTDIR <output> 
        
        --BASE Run basename for organization
        --OUTDIR Output directory
        --CHAIN  Liftover chain file 
        --BED Input bed file 
        --REF_FASTA reference fasta genome for bedtools getfasta extraction
    
    """.stripIndent()
}