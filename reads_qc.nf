#!/usr/bin/env nextflow

//NXF_ANSI_LOG=false
//script parameters
params.reads	=  "/home/drewx/Documents/ez-pool-seq/RawReads"
params.pattern 	=  "*_R{1,2}_001.fastq.gz"
params.hcpu	=  4
reads_pattern 	=  params.reads + "/*" + params.pattern

Channel.fromFilePairs(reads_pattern)
       .ifEmpty{ exit 1, "params.reads empty no reads found" }
       .into{ raw_reads; raw_reads_trimgalore }

raw_reads.map{ it  -> [ it[1][0], it[1][1]] }
         .set{ raw_reads_FastQC }


process FastQC{

    echo true
    publishDir "$PWD/FastQC" 
    cpus params.hcpu 
    input:
	file reads from raw_reads_FastQC.collect()

	
    output:
	set  file("*fastqc.html"),  file("multiqc_*") into FastQC_results

    script:


"""

    fastqc \
        -format fastq \
        -threads $params.hcpu \
        $reads
	
    multiqc .
  
"""
     	
}


