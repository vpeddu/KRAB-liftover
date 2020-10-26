OUTDIR = "${params.OUTDIR}"
BEDFILE = file("${params.BED}")
CHAINFILE = file ("${params.CHAIN}")
REF_FASTA = file("${params.REFERENCE}")
BASE = "${params.base}"


process liftOver {  
    input: 
        file BEDFILE 
        file CHAINFILE
        val BASE
    container 'quay.io/biocontainers/ucsc-liftover:357--h446ed27_4'   
    publishDir "${OUTDIR}/liftOver", mode: 'copy'

    output:
        file "${BASE}.transferred.bed" into liftoverOutcH
        file "${BASE}.unmapped.bed"
    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    liftOver ${BEDFILE} ${CHAINFILE} ${BASE}.transferred.bed ${BASE}.unmapped.bed

    """
}

process getFasta {  
    input: 
        file BEDFILE from liftoverOutcH
        file REF_FASTA
        val BASE
    container 'quay.io/biocontainers/bedtools:2.26.0gx--he513fc3_4'   
    publishDir "${OUTDIR}/getFasta", mode: 'copy'

    output:
        file "${BASE}.extracted.fasta" into toSmashCh

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    bedtools getfasta -s -name -fi ${REF_FASTA} -bed ${BEDFILE} > ${BASE}.extracted.fasta


    """
}

process smashFasta {  
    input: 
        file EXONFASTA from toSmashCh
        val BASE

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