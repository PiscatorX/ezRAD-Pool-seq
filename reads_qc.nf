#!/usr/bin/env nextflow

//NXF_ANSI_LOG=false
//script parameters
params.reads	=  "/media/drewx/1D995C4565A1FA44/Users/A0047731/My Documents/dx.tmp/Pool-seq"
params.pattern 	=  "R{1,2}_001.fastq.gz"

reads_pattern 	= file(params.reads + "/*" + params.pattern)
raw_reads     	= Channel.fromFilePairs(reads_pattern)


process FastQC{
    echo true
    input:
	set read1, read2 from raw_reads

    output:
	
	
"""

echo $read1 

"""
     	
}