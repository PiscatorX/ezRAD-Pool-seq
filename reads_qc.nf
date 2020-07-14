#!/usr/bin/env nextflow

NXF_ANSI_LOG=false
//script parameters
params.reads	=  "/home/drewx/Documents/ez-pool-seq/RawReads2"
params.pattern 	=  "*_R{1,2}_001.fastq.gz"

//params.pattern 	=  "*_{1,2}.fastq"
params.output   =  "$PWD/ezRAD_Test"
params.hcpu	=  4
params.bwa_ref  =  "/opt/DB_REF/Genomes/Zostera_marina/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna"


reads_pattern 	=  params.reads + "/*" + params.pattern

Channel.fromFilePairs(reads_pattern)
       .ifEmpty{ exit 1, "params.reads empty no reads found" }
       .into{raw_reads_trimgalore; raw_reads_bwa }
       //raw_reads;
       

//raw_reads.map{ it  -> [ it[1][0], it[1][1]] }
//         .set{ raw_reads_FastQC }


// process FastQC{

//     echo true
//     publishDir "$output/FastQC" 
//     cpus params.hcpu 
//     input:
// 	file reads from raw_reads_FastQC.collect()

	
//     output:
// 	set  file("*fastqc.html"),  file("multiqc_*") into FastQC_results

//     script:


// """

//     fastqc \
//         -format fastq \
//         -threads ${params.hcpu} \
//         $reads
	
//     multiqc .
  
// """
     	
// }




// process Trim_Galore{

//     echo true
//     publishDir params.output + "/TrimGalore" 
//     cpus params.hcpu 
//     input:
//          set val(sample),  file(reads) from raw_reads_trimgalore

//     output:
// 	 file("*val_{1,2}.fq.gz") into TrimGalore
// 	 file("TrimQC") into Trim_galore_FastQC
	 
//     script:


// """
//     mkdir TrimQC
    
//     trim_galore \
//         --paired \
//         --cores $params.hcpu \
//         --quality 20 \
//         --length  30 \
//         --fastqc \
//         $reads

//     multiqc  .
    
//     mv *fastqc.html  *report*   TrimQC
	      
// """
     	
// }




process  BWA_MEM{

       echo true
       publishDir params.output + "/BWA-MEM/$sample" 
       cpus params.hcpu 
       input:
            set val(sample),  file(reads) from raw_reads_bwa

      output:
           file "${sample}.*"  into SAM_files
	   file "Plots"  into stat_plots

"""
   
   bwa mem \
       $params.bwa_ref \
       $reads \
       -t ${params.hcpu} \
       -o ${sample}.sam \
       -T 20

   samtools \
       view \
       -t ${params.hcpu} \
       ${sample}.sam \
       -o ${sample}.bam 
            
   samtools \
      sort \
      -t ${params.hcpu} \
      ${sample}.bam \
      -o ${sample}.sorted.bam

   samtools \
      index \
      ${sample}.sorted.bam

   samtools \
      stats \
      ${sample}.sorted.bam > \
      ${sample}.sorted.stats

  samtools \
      idxstats \
      --threads ${params.hcpu} \
      ${sample}.sorted.bam > \
      ${sample}.sorted.idxstats

  plot-bamstats \
      ${sample}.sorted.stats \
      -p Plots/${sample}
      
  multiqc  . \
      --outdir \
      Plots


"""
//idxstats sequence name, sequence length, # mapped read-segments and # unmapped read-segments.
}

//TO DO 
// Mapped reads were subsampled to median coverage in SAMTools using the view command with the ‘-s’ flag.
//Struggling to figure how this out
//Median coverage of all reads combined?
//How was this done?
//Custom script
//Will look at this when I have SNP counts