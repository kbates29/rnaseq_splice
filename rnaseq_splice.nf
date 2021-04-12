#!/vcu_gpfs2/home/batesk2/bin/nextflow

/*
 * Starting with a raw fastq performs adapter trimming and QC with Trim Galore, 
 * aligns to reference genome with STAR and obtains gene counts with featurecounts
 */


/*
 * Defines parameters
 */

params.reads = null
params.genome_dir = null
params.gtf = null
params.out_dir = null
params.bed = null
params.group1 = null
params.group2 = null
params.name1 = null
params.name2 = null

/*
 * Create channel for reads
 */
 Channel
    .fromFilePairs(params.reads)
    .ifEmpty {error "Cannot find reads matching: ${params.reads}"}
    .into {reads_fast_qc; reads_trim}


/*
 * Run fastqc prior to trimming to compare afterwords
 */
process fast_qc{
    executor 'sge'
    clusterOptions '-S /bin/bash -V'
    tag "$sampleId"
    publishDir "${params.out_dir}/fastqc", mode: 'copy', 
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$sampleId/$filename" : "$filename"}

    input: 
        tuple sampleId, file(reads) from reads_fast_qc
    
    output:
        file "*_fastqc.{zip,html}" into fastqc_results

    script:
        """
	    fastqc -q $reads
        """
}

/*
 *  Run Trim Galore which trims adapters with cutadapt and then reruns fastqc
 */

process trim_galore{
    executor "sge"
    clusterOptions="-S /bin/bash -V -pe smp 8"
    tag "$sampleId"
    publishDir "${params.out_dir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$sampleId/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$sampleId/$filename"
        }

    input:
        set sampleId, file(reads) from reads_trim

    output:
        file "*.gz" into trimmed_reads
        file "*trimming_report.txt" into trimgalore_results
        file "*_fastqc.{zip,html}" into trimgalore_fastqc
    
    script:
        """
  	    trim_galore --paired -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --cores 4 --fastqc $reads
        """
}

/*
 * Run star and index bam files
 */
 process star{
    executor "sge"
    clusterOptions="-S /bin/bash -V -pe smp 8"
    tag "$sampleId"
    publishDir "${params.out_dir}/star", mode: 'copy',
        saveAs: { filename ->
            if (filename.indexOf("Log") >0 ) "logs/$filename"
            else if (filename.indexOf("STARgenome") >0) "STARgenome/$filename"
            else if (filename.indexOf("STARpass1") >0) "STARpass1/$filename"
            else ("$filename")
        }
       
    input:
        file(reads) from trimmed_reads

    output:
        file "*.bam" into star_aligned, junc_ann, junc_sat, read_dist,create_bam, bam_big
        file "*SJ.out.tab" into splice_junction
        file "*"
        file "*Log*.out" into star_log
        file "*STARgenome"
        file "*STARpass1"
        file "*.bai" into index_bam, ann_bai, sat_bai, dist_bai, coverage_bai
        

    script:
        sampleId = reads[0].toString() - ~/(_1)?(_2)?(\.txt)?(\.gz)?(_trimmed)?(_val_1)?(_val_2)?(\.fq)?(\.fastq)?(\.txt)?(\.gz)?$/
        sampleId = sampleId + "_"
        
        """        
        STAR --runMode alignReads --runThreadN 8 \
        --genomeDir ${params.genome_dir} \
        --readFilesIn $reads \
        --sjdbGTFfile ${params.genome_dir}/${params.gtf} \
        --sjdbOverhang 100 \
        --outFileNamePrefix $sampleId \
        --outSAMtype BAM SortedByCoordinate \
        --twopassMode Basic --readFilesCommand zcat

        samtools index ${sampleId}Aligned.sortedByCoord.out.bam
        """
 }

/*
 * Run Junction Annotation
 */
 process junction_annotation {
    executor "sge"
    clusterOptions="-S /bin/bash -V"
    tag "$sampleId"
    publishDir "${params.out_dir}/rseqc" , mode: 'copy',
          saveAs: {filename ->
                   if (filename.indexOf("junction_plot.r") > 0)     "junction_annotation/rscripts/$filename"
              else if (filename.indexOf("junction.xls") > 0)        "junction_annotation/data/$filename"
              else if (filename.indexOf("splice_events.pdf") > 0)   "junction_annotation/events/$filename"
              else if (filename.indexOf("splice_junction.pdf") > 0) "junction_annotation/junctions/$filename"
              else filename
          }
   
   input:
        file bam from junc_ann
        file bai from ann_bai
   
   output:
        file "*.{txt,pdf,r,xls}" into ann_results
   
   script:
    sampleId = bam[0].toString() - '_Aligned.sortedByCoord.out.bam'
   """
    junction_annotation.py -i $bam -o $sampleId -r ${params.genome_dir}/${params.bed} 2> ${sampleId}.junction_annotation_log.txt
   """
 }


 /*
  * Run junction_saturation
  */
  process junction_saturation {
    executor "sge"
    clusterOptions="-S /bin/bash"
    tag "$sampleId"
    publishDir "${params.out_dir}/rseqc" , mode: 'copy',
          saveAs: {filename ->
                 if (filename.indexOf("junctionSaturation_plot.pdf") > 0)    "junction_saturation/$filename"
              else if (filename.indexOf("junctionSaturation_plot.r") > 0)    "junction_saturation/rscripts/$filename"
              else filename
          }

    input:
        file bam from junc_sat
        file bai from sat_bai
    
    output:
        file "*.{txt,pdf,r,xls}" into sat_results

    script:
    sampleId = bam[0].toString() - '_Aligned.sortedByCoord.out.bam'
    """
    junction_saturation.py -i $bam -r ${params.genome_dir}/${params.bed} -o $sampleId
    """
  }

/*
 * Run read_disctribution
 */
 process read_distribution {
    executor "sge"
    clusterOptions="-S /bin/bash"
    tag "$sampleId"
    publishDir "${params.out_dir}/rseqc" , mode: 'copy',
          saveAs: {filename ->
                if (filename.indexOf("read_distribution.txt") > 0) "read_distribution/$filename"
                else filename

          }
   input:
        file bam from read_dist
        file bai from dist_bai
   
   output:
        file "*.txt" into dist_results
   script:
   sampleId = bam[0].toString() - '_Aligned.sortedByCoord.out.bam'
   """
   read_distribution.py  -i $bam -r ${params.genome_dir}/${params.bed} > ${sampleId}.read_distribution.txt
   """
 }



/*
 * Create bigwig file for viewing
 */
 process bamcoverage{
    executor "sge"
    clusterOptions="-S /bin/bash -V -pe smp 8"
    tag "$sampleId"
    publishDir "${params.out_dir}/bigwig", mode: 'copy'

    input: 
        file bam from bam_big
        file bai from coverage_bai
 
    output:
        file "*.bw" into bigwig 
 
    script:
    sampleId = bam[0].toString() - '_R1_Aligned.sortedByCoord.out.bam'
    """
    bamCoverage -b $bam -o ${sampleId}.bw -p 8 -bs 10
    """
 }

 /*
 *  Run featurecounts 
 */
 process featurecounts{
    executor "sge"
    clusterOptions="-S /bin/bash -V -pe smp 8"
    tag "$sampleId"
    publishDir "${params.out_dir}/counts", mode: 'copy'
 
    input:
        file bam from star_aligned

    output:
        file "*counts.txt" into sample_counts
        file "*.summary" into count_log

    script:
        sampleId = bam[0].toString() - 'Aligned.sortedByCoord.out.bam'
        """
        featureCounts -p -T 8 -t exon -g gene_id -s 2 -a ${params.genome_dir}/${params.gtf}\
        -o ${sampleId}counts.txt $bam
        """
 }


/*
 * Merge counts into one file 
 */
 process mergecounts{
    executor "sge"
    clusterOptions="-S /bin/bash -V"
    publishDir "${params.out_dir}/counts", mode: 'copy'   
   
    input:
        file input_files from sample_counts.collect()
   
    output:
        file 'merged_counts.txt'

    script:
        gene_id = "<(tail -n +2 ${input_files[0]} | cut -f1)"
        counts = input_files.collect{filename ->
            "<(tail -n +2 ${filename} | cut -f7)"}.join("\t")
        """
        paste $gene_id \t $counts > merged_counts.txt
        """
 }
 

/*
 * MultiQC
 */
 process multiqc {
     executor "sge"
     clusterOptions = "-S /bin/bash -V"
     publishDir "${params.out_dir}/multiqc", mode: 'copy'

    input:
        file fastqc from fastqc_results.collect()
        file star from star_log.collect()
        file sat from sat_results.collect()
        file ann from ann_results.collect()
        file dist from dist_results.collect()
        file counts from count_log.collect()
        file trimming from trimgalore_results.collect()
        file trimqc from trimgalore_fastqc.collect()

    output:
        file "*.html" into multiqc_report
        file "*" into mutli_data
    
    script:
        """
        multiqc .
        """

 }


/*
 * Sort bams into groups
 */
 process bam_group {
    executor "sge"
    clusterOptions="-S /bin/bash -V"
    publishDir "${params.out_dir}/star", mode:'copy'

    input: 
        file bams from create_bam.collect()
    
    output:
        file "*.txt" into rmats_groups
    
    script:
        test = bams.join(",")
        """
        python /vcu_gpfs2/home/batesk2/python_scripts/rmats_groups.py -d ${params.out_dir}/star/ -b $test -g1 ${params.group1} \
        -g2 ${params.group2} -n1 ${params.name1} -n2 ${params.name2}
        """



 }


  /*
  * rMATS prep step per file
  */
  process rmats {
    executor "sge"
    clusterOptions="-S /bin/bash -V -pe smp 8"
    publishDir "${params.out_dir}/rmats", mode:'copy'
    
    input:
        file (bams) from rmats_groups.collect()
  
    output:
        file "*.log" into rmats_log

    script:
    bam_1 = bams[0]
    bam_2 = bams[1]
    """
        python /vcu_gpfs2/home/batesk2/bin/rmats-turbo/rmats.py --b1 $bam_1 --b2 $bam_2  \
        --gtf ${params.genome_dir}/${params.gtf} -t paired  --readLength 75 --nthread 8 \
        --od ${params.out_dir}/rmats --tmp ${params.out_dir}/rmats_tmp > rmats.log
    """
  }

