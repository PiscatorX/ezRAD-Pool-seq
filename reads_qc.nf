#!/usr/bin/env nextflow

NXF_ANSI_LOG=false
//script parameters
params.reads	=  "/home/drewx/Documents/ez-pool-seq/RawReads2"
params.pattern 	=  "*_R{1,2}_001.fastq.gz"

//params.pattern 	=  "*_{1,2}.fastq"
params.output   =  "$PWD/Clinid_/"
params.hcpu	=  4
params.bwa_ref  =  "/opt/DB_REF/Genomes/Zostera_marina/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna"
params.faidx_ref  =  "/opt/DB_REF/Genomes/Zostera_marina/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna.fai"


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

       publishDir params.output + "/bwa_mem/sam/" 
       cpus params.hcpu 
       input:
            set val(sample),  file(reads) from raw_reads_bwa

      output:
           set val(sample), file("${sample}.*")  into SAM_files

"""

     bwa mem \
	 $params.bwa_ref \
	 $reads \
	 -t ${params.hcpu} \
	 -o ${sample}.sam \
	 -T 20
"""

}



process  samtools_index{

       publishDir params.output + "/BAM/" 
       cpus params.hcpu 
       input:
            set val(sample),  file(sam_file) from SAM_files

      output:
           set val(sample), file("${sample}_sorted.bam") into (bam_stats, bam_RG)
	   set val(sample), file("${sample}_sorted.bam.bai") into (index_files)
	   file("${sample}_sorted.bam.bai") into index_files2
	   
"""
     samtools \
	 view \
	 -@ ${params.hcpu} \
	 ${sample}.sam \
	 -o ${sample}.bam 

     samtools \
	sort \
	-@ ${params.hcpu} \
	${sample}.bam \
	-o ${sample}_sorted.bam

     samtools \
	index \
        -@ ${params.hcpu} \
	${sample}_sorted.bam

"""
//http://www.htslib.org/doc/samtools.html

}


process  bam_stats{

       
       publishDir params.output + "/BAM_stats" 
       cpus params.hcpu 
       input:
            set val(sample),  file(bam_file) from bam_stats
	    set val(sample), file(index) from index_files

      output:
           file "${sample}_sorted*"  into stats
	   file "plots"  into stat_plots

"""

     samtools \
    	stats \
    	${sample}_sorted.bam > \
    	${sample}_sorted.stats

    samtools \
    	idxstats \
    	--threads ${params.hcpu} \
    	${sample}_sorted.bam > \
    	${sample}_sorted.idxstats

    plot-bamstats \
    	${sample}_sorted.stats \
    	-p plots/${sample}

    multiqc  . \
    	--outdir \
    	plots
    

"""
//idxstats sequence name, sequence length, # mapped read-segments and # unmapped read-segments.
}
//to do 
// mapped reads were subsampled to median coverage in samtools using the view command with the ‘-s’ flag.
//struggling to figure how this out
//median coverage of all reads combined?
//how was this done?
//custom script
//will look at this when i have snp counts



process  AddReadGroups{

     
     publishDir params.output + "/BAM/"
     
     input:
           set val(sample), file(bam_file) from bam_RG

     output:
           file("${sample}_rg.bam") into (bam_mpileup1, bam_mpileup2)

"""

   $picard AddOrReplaceReadGroups \
       I=${bam_file} \
       O=${sample}_rg.bam \
       RGID=${sample} \
       RGLB=lib1 \
       RGPL=ILLUMINA \
       RGPU=unit1 \
       RGSM=${sample}


"""

}


Channel.fromPath("/home/drewx/Documents/ez-pool-seq/ClinidX/BAM/*sorted.bam").into{bam_mpileup1; bam_mpileup2}

process mpileup{
  
    publishDir params.output + "/mpileup"
    input:
         file bam_files from bam_mpileup1.collect()
	
    output:
        file("bam_files")
	file("samtools.pileup") into samtools_mpileup
	    
"""

    ls *.bam > bam_files

    samtools \
        mpileup \
        --fasta-ref ${params.bwa_ref} \
        --min-BQ 20 \
        --max-depth 1000 \
        --bam-list bam_files \
        -o samtools.pileup
      
 
"""
    
//http://www.htslib.org/doc/samtools-mpileup.html    

}
//TO DO
//Merge BAM files for igv



process bcftools{

    publishDir params.output + "/bcftools"
    input:
        file bam_files from bam_mpileup2.collect()
        

    output:
        file("bam_files")
	file("bcftools_calls.vcf") into bcftools_mpileup
	    
"""

    ls *.bam > bam_files
      
    bcftools \
       mpileup \
       --bam-list bam_files \
       --max-depth 1000 \
       --min-BQ 20 \
       -f ${params.bwa_ref} \
       alignments.bam | bcftools \
       call \
       -mv \
       -o bcftools_calls.vcf

  
"""
//https://samtools.github.io/bcftools/howtos/variant-calling.html

}




process pileup2SNP{


    echo true
    //errorStrategy 'ignore'
    publishDir params.output + "/popooltn2"
    
    input:
         file pileup from samtools_mpileup
	 
    output:
        set file("mpileup.sync"), file("mpileup_SNP*") into pileup_snps


"""

    mpileup2sync.pl \
        --fastq-type sanger \
        --input ${pileup} \
        --output mpileup.sync
   
    snp-frequency-diff.pl \
        --input mpileup.sync \
        --output-prefix mpileup_SNP \
        --min-count 4 \
        --min-coverage 10 \
        --max-coverage 500

"""
//TO DO
//Must check the file order in the file    
//mpileup2sync="java -ea -Xmx7g -jar /opt/popoolation2_1201/mpileup2sync.jar --threads ${params.hcpu}"
//https://sourceforge.net/p/popoolation2/wiki/Tutorial/
}










