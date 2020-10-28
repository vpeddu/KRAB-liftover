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

    cpus 1
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
        tuple(val(BASE), file("${BASE}.smashed.fasta")) into translateCh

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    Rscript --vanilla ${baseDir}/bin/smash_fasta.r ${EXONFASTA}

    mv smashed.fasta ${BASE}.smashed.fasta

    """
}

process faTrans {  
    input: 
        tuple(val(BASE), file(SMASHED_FASTA)) from translateCh
        //val BASE

    container 'genomicpariscentre/kentutils'   
    publishDir "${OUTDIR}/faTrans", mode: 'copy'

    output:
        tuple(val(BASE), file("${BASE}.translated.fasta")) into alignCh

    cpus 1
    memory 2.Gb 

    script:
    """
    #!/bin/bash

    faTrans ${SMASHED_FASTA} ${BASE}.translated.fasta


    """
}


process MAFFT {  
    input: 
        tuple(val(BASE), file(TRANSLATED_FASTA)) from alignCh
        //val BASE

    container 'staphb/mafft'   
    publishDir "${OUTDIR}/MAFFT", mode: 'copy'

    output:
        tuple(val(BASE),file("${BASE}.aligned.fasta")) into treeCh

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    mafft --thread ${task.cpus} ${TRANSLATED_FASTA} > ${BASE}.aligned.fasta


    """
}

process fastTree {  
    input: 
        tuple(val(BASE), file(ALIGNED_FASTA)) from treeCh
        //val BASE

    container 'staphb/fasttree'   
    publishDir "${OUTDIR}/fastTree", mode: 'copy'

    output:
        tuple(val(BASE),file("${BASE}.tree.newick")) 
        file "${BASE}.tree.newick" into compareTreesCh

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    export OMP_NUM_THREADS=${task.cpus}

    FastTree ${ALIGNED_FASTA} > ${BASE}.tree.newick


    """
}

//compareTreesCh.toList().view()

process compareTrees {  
    input: 
        file treefiles from compareTreesCh.collect()
        //val BASE

    container 'quay.io/biocontainers/r-phytools:0.6_99--r40h6115d3f_1'   
    publishDir "${OUTDIR}", mode: 'copy'

    output:
        file "comparison_info.csv"

    cpus 2
    memory 4.Gb 

    script:
    """
    #!/bin/bash

    Rscript --vanilla  ${baseDir}/bin/compare_trees.r 


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