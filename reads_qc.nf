#!/usr/bin/env nextflow

NXF_ANSI_LOG=false
//script parameters
params.reads  =  "/home/drewx/Documents/ez-pool-seq/Z_marina/"
params.pattern 	=  "*{1,2}.fastq.gz"
params.target_coverage = 30
params.pool_size = params.target_coverage
params.max_coverage = 500
params.min_coverage = 5
params.mia_count    = 2
params.min_qual = 20
params.output   =  "$PWD/Z_capensis/"
params.hcpu	=  4

params.bwa_ref  =  "/opt/DB_REF/Genomes/Zostera_marina/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna"
// params.faidx_ref  =  "/opt/DB_REF/Genomes/Zostera_marina/GCA_001185155.1_Zosma_marina.v.2.1_genomic.fna.fai"
// params.gtf =  "/opt/DB_REF/Genomes/Zostera_marina/GCA_001185155.1_Zosma_marina.v.2.1_genomic.gtf2"


reads_pattern 	=  params.reads + "/*" + params.pattern

Channel.fromFilePairs(reads_pattern)
       .ifEmpty{ exit 1, "params.reads empty no reads found" }
       .into{raw_reads; raw_reads_trimgalore}
       

raw_reads.map{ it  -> [ it[1][0], it[1][1]] }
         .set{ raw_reads_FastQC }



process FastQC{

    publishDir  params.output + "/FastQC", mode: 'move' 
    cpus params.hcpu 
    input:
	file reads from raw_reads_FastQC.collect()

	
    output:
	set  file("*fastqc.html"),  file("multiqc_*") into FastQC_results

    script:


"""

    fastqc \
        -format fastq \
        -threads ${params.hcpu} \
        $reads
	
    multiqc .
  
"""
     	
}




process Trim_Galore{

    publishDir params.output + "/TrimGalore", mode: 'copy'
    cpus params.hcpu 
    input:
         set val(sample),  file(reads) from raw_reads_trimgalore

    output:
	 set val(sample), file("*val_{1,2}.fq.gz") into  reads_bwa
	 file("TrimQC") into Trim_galore_FastQC
	 
script:
"""
    mkdir TrimQC
    
    trim_galore \
        --paired \
        --cores $params.hcpu \
        --quality 20 \
        --length  30 \
        --fastqc \
        $reads

    multiqc  .
    
    mv *fastqc.html  *report*   TrimQC
	      
"""
     	
}



process  BWA_MEM{

       
     publishDir params.output + "/bwa_mem/sam/", mode: 'copy'
       cpus params.hcpu 
       input:
         set val(sample),  file(reads) from reads_bwa

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
//TO DO
//Add Read groups




process  samtools_index{

     publishDir params.output + "/BAM/" , mode: 'copy'
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

       
    publishDir params.output + "/BAM_stats", mode: 'copy'
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

//TO DO 
// mapped reads were subsampled to median coverage in samtools using the view command with the ‘-s’ flag.
//struggling to figure how this out
//median coverage of all reads combined?
//how was this done?
//custom script
//will look at this when i have snp counts



process  AddReadGroups{

     
    publishDir params.output + "/BAM/", mode: 'copy'
     
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



process mpileup{
  
    publishDir params.output + "/mpileup", mode: 'copy'
    input:
        file bam_files from bam_mpileup1.collect()
	
    output:
        file("bam_files") into sample_list
	file("samtools.pileup") into (samtools_mpileup1,
                                      samtools_mpileup2,
				      samtools_mpileup3,
				      samtools_mpileup4)
				      
	    
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

    publishDir params.output + "/bcftools", mode: 'move'
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


    publishDir params.output + "/popooltn2", mode: 'copy'
    
    input:
         file pileup from samtools_mpileup1
	 
    output:
        file("mpileup.sync") into mpileup_sync
        file("mpileup_SNP_rc") into (mpileup_snps_rc1, mpileup_snps_rc2) 
	   file("mpileup_SNP_pwc") into mpileup_snps_pwc


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



process rc2region{


    echo true
    publishDir params.output + "/popooltn2", mode: 'copy'
    input:
         file mpileup_snps_rc1
	 

    output:
         file("regions.txt") into (genome_region1, genome_region2)

"""

    get_region.py \
        --target_coverage ${params.target_coverage} \
        ${mpileup_snps_rc1} > regions.txt


"""

}



process splitRC{


    echo true
    publishDir params.output + "/sample_SNP", mode: 'copy'
    input:
         file mpileup_snps_rc2
	 file sample_list

    output:
         file("*.snps") into snp_files

"""

    splitRC.py \
       ${mpileup_snps_rc2} \
       -s ${sample_list}

"""

}



    
process countSNPs{

    echo true
    publishDir params.output + "/popooltn2", mode: 'move'
    input:
        file sample_snp  from snp_files.collect()

    output:
         file("SNP.counts")
   
script:
"""

    grep -c "snp\$" ${sample_snp} | tee  SNP.counts  


"""

}



process sync2GenePop{

    publishDir params.output + "/popooltn2", mode: 'copy'
    input:
         each region from genome_region1.splitText()
	 file sync_file from mpileup_sync
	 
    output:
        file("${region_id}.genepop") into genepop_merge
	  
	     

script:
(region_id,coord) = region.tokenize(":")


"""

    subsample_sync2GenePop.pl \
         --input ${sync_file} \
         --output ${region_id}.genepop \
         --target-coverage ${params.target_coverage} \
         --max-coverage ${params.max_coverage} \
         --method fraction \
          --diploid \
         --region ${region}
                 
"""

}



process Genepop_merge{

    //echo true
    publishDir params.output + "/popooltn2", mode: 'copy'
    input:
         file genepop from genepop_merge.collect()
	 
    output:
        file("pool.Genepop") into pool
	file("Genepop.ref") into ref
	  	     

script:
"""
    merge_genepop.py *.genepop     
                 
"""

}
 


process TajimaPi{

    publishDir params.output + "/diversityMetrics", mode: 'copy'
    input:
         file samtools_mpileup2

    output:
         file("${samtools_mpileup2.baseName}.*")

	 
         
"""

   Variance-sliding.pl \
       --input ${samtools_mpileup2} \
       --output ${samtools_mpileup2.baseName}.pi \
       --snp-output ${samtools_mpileup2.baseName}.pi.snps \
       --fastq-type sanger \
       --measure pi \
       --window-size 100 \
       --step-size 100 \
       --pool-size ${params.pool_size} \
       --min-qual ${params.min_qual} \
       --min-count ${params.mia_count} \
       --min-coverage ${params.min_coverage} \
       --max-coverage ${params.max_coverage} 

           
"""


}




process TajimaD{

    publishDir params.output + "/diversityMetrics", mode: 'copy'
    input:
         file samtools_mpileup4

    output:
         file("${samtools_mpileup4.baseName}.*")

	 
         
"""

   Variance-sliding.pl \
       --input ${samtools_mpileup4} \
       --output ${samtools_mpileup4.baseName}.D \
       --snp-output ${samtools_mpileup4.baseName}.D.snps \
       --fastq-type sanger \
       --measure D \
       --window-size 100 \
       --step-size 100 \
       --pool-size ${params.pool_size} \
       --min-qual ${params.min_qual} \
       --min-count ${params.mia_count} \
       --min-coverage ${params.min_coverage} \
       --max-coverage ${params.max_coverage} 
  
           
"""
}




process WattersonTheta{

    publishDir params.output + "/diversityMetrics", mode: 'copy'
    input:
         file samtools_mpileup4

    output:
         file("${samtools_mpileup4.baseName}.*")
         
"""

   Variance-sliding.pl \
       --input ${samtools_mpileup4} \
       --output ${samtools_mpileup4.baseName}.theta \
       --snp-output ${samtools_mpileup4.baseName}.theta.snps \
       --fastq-type sanger \
       --measure theta \
       --window-size 100 \
       --step-size 100 \
       --pool-size ${params.pool_size} \
       --min-qual ${params.min_qual} \
       --min-count ${params.mia_count} \
       --min-coverage ${params.min_coverage} \
       --max-coverage ${params.max_coverage} 
           
"""


}
