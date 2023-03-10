/* 
 * pipeline input parameters 
 */

params.genome = "/Zulu/bnolan/Indexes/bwaIndex/hg38.fa" 
params.samplesheet = "${baseDir}/Samplesheets/samples_test.csv"
params.outdir = "${baseDir}/results"
params.index = "/Zulu/bnolan/Indexes/bwaIndex/"
params.threads = "4"
params.aligner = "bwa-meth"


log.info """\
         ===================================
         B S S E Q - N F   P I P E L I N E    
         ===================================
         genome       : ${params.genome}
         outdir       : ${params.outdir}
         samplesheet  : ${params.samplesheet}
         threads      : ${params.threads}
         aligner      : ${params.aligner}

         ${params.aligner}
         ------------------------------------
         index        : ${params.index}

         """
         .stripIndent()


// Parse samplesheet and create reads channel            
Channel
        .from ( file(params.samplesheet) )
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]] }
        .set { read_pairs_ch }

        
index_ch = Channel.value(file(params.index))

genome_ch = Channel.value(file(params.genome)) 


// // Concatenate reads of same sample [e.g. additional lanes, resequencing of old samples]
//nfcore
read_pairs_ch
        .groupTuple()
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }  


process CAT_FASTQ {
    //nf-core
    tag "$sample_id"
    publishDir "${params.outdir}/fastqMerged", mode: 'copy'

    input:
    tuple val(sample_id), path(reads, stageAs: "input*/*") from ch_fastq.multiple

    output:
    tuple val(sample_id), path("*.merged.fastq.gz") into cat_out_ch


    script:
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]

     if (readList.size >= 2) {
        def read1 = []
        def read2 = []
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }

        """
        cat ${read1.join(' ')} > ${sample_id}_1.merged.fastq.gz
        cat ${read2.join(' ')} > ${sample_id}_2.merged.fastq.gz
        """  
     }
}


//mix merged reads with singles
cat_out_ch
        .mix(ch_fastq.single)
        .into { cat_merged_ch; cat_merged_ch2 }



//  Run fastQC to check quality of reads files

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", pattern:"{*.html,fastqc_${sample_id}_logs}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from cat_merged_ch

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}


// Trimming reads with Trim Galore

process trimming {
    tag "$key"
    publishDir "${params.outdir}/trimmed", pattern: "*.fq.gz", mode: 'copy'

    input: 
    tuple val(key), path(reads) from cat_merged_ch2

    output:
    //tuple path(fq_1_paired), path(fq_2_paired) into ch_out_trimmomatic
    tuple val(key), path("*.fq.gz") into ch_out_trimmomatic, ch_out_trimmomatic2

    script:

    """
    trim_galore \\
                --paired \\
                ${reads[0]} \\
                ${reads[1]} \\
                --clip_R1 5 \\
                --clip_R2 5 \\
                --three_prime_clip_R1 5 \\
                --three_prime_clip_R2 5 \\
                --basename $key \\
                --cores $params.threads
    """
}

//  Run fastQC to check quality of reads files

process fastqc_trimmed {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_trimmed", pattern:"{*.html,fastqc_${sample_id}_trimmed_logs}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from ch_out_trimmomatic

    output:
    path("fastqc_${sample_id}_trimmed_logs") into fastqc_trimmed_ch

    script:
    """
    mkdir fastqc_${sample_id}_trimmed_logs
    fastqc -o fastqc_${sample_id}_trimmed_logs -f fastq -q ${reads}
    """  
}



// Align reads to index with Bismark or bwa-meth

process align {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/alignment/$key/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(reads) from ch_out_trimmomatic2
    path(index) from index_ch   // index_ch

    output:
    tuple val(key), path("*bam") into bam_ch, bam_ch2
    path("*report.txt"), optional: true into align_report_ch // Not used if params.aligner == bwa-meth

    script:

    if( params.aligner == 'bismark' ) 
        """
        bismark --genome $index \\
                -1 ${reads[0]} \\
                -2 ${reads[1]} \\
                --bam \\
                --multicore $params.threads
        """
    else if( params.aligner == 'bwa-meth')
        """
        INDEX=`find -L ${index} -name "*.bwameth.c2t" | sed 's/.bwameth.c2t//'`
        # Modify the timestamps so that bwameth doesn't complain about building the index
        # courtesy of nf-core
        #touch -c -- *

        bwameth.py \\
                -t $params.threads \\
                --reference \$INDEX \\
                $reads \\
                | samtools view -@ $params.threads -bhS -o ${key}.bam -
        """   
}


// Samtools stats for summary statistics on bwa-meth alignment
process samtools_stat_flagstat {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/samtools/stats/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_ch

    when:
    params.aligner == 'bwa-meth'

    output:
    path('*stats') into samtools_stats_ch
    path('*flagstat') into samtools_flag_ch


    script:
    """
    samtools stats --threads ${params.threads} \\
                   $bam \\
                   > ${key}.stats

    samtools flagstat \\
                        --threads ${params.threads} \\
                        $bam \\
                        > ${key}.flagstat

    """
}


//deduplication
//sort
process samtools_sort {
        tag "$key"
        publishDir "${params.outdir}/${params.aligner}/samtools/", pattern:'*', mode: 'copy'

        input:
        tuple val(key), path(bam) from bam_ch2
        
        output:
        tuple val(key), path('*sorted.bam') into bam_sorted_ch
        tuple val(key), path('*sortn.bam') into bam_sortn_ch
        
        script:
        """
        samtools sort -n $bam > ${key}.sortn.bam
        samtools sort $bam > ${key}.sorted.bam
        """
    }


// Deduplicate the bam files with bismark
process deduplication_bismark {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/deduplicate/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sortn_ch

    when:
    params.aligner == 'bismark'

    output:
    tuple val(key), path("*.deduplicated.bam") into dedup_bam_bismark_ch
    path("*.deduplication_report.txt") into dedup_report_bismark_ch

    script:
    """
    deduplicate_bismark -p --bam $bam
    """
}



// Deduplicate bam files with picard [bwa-meth]
process deduplication_picard{
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/picard/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sorted_ch 
    path(genome) from genome_ch

    when:
    params.aligner == 'bwa-meth'

    output:
    tuple val(key), path("*.deduplicated.bam") into dedup_bam_picard_ch
    path("*.MarkDuplicates.metrics.txt") into dedup_report_picard_ch

    script:
    """
    picard \\
        -Xmx3g \\
        MarkDuplicates \\
        --INPUT $bam \\
        --OUTPUT ${key}.deduplicated.bam \\
        --REFERENCE_SEQUENCE $genome \\
        --REMOVE_DUPLICATES true \\
        --METRICS_FILE ${key}.MarkDuplicates.metrics.txt
    """
}


if ( params.aligner == 'bismark'){
    dedup_bam_bismark_ch
                    .into{ dedup_bam_ch; dedup_bam_ch2 }
    dedup_report_ch = dedup_report_bismark_ch
}else{
    dedup_bam_picard_ch
                    .into{ dedup_bam_ch; dedup_bam_ch2 }
    dedup_report_ch = dedup_report_picard_ch
} 


//Combine replicates: split by '_', and group all samples
dedup_bam_ch
    .map { group_rep, bam ->
                        def(group) = group_rep.split("_")  
                        tuple( group, bam )
                        }
    .groupTuple()
    .set { dedup_groupsplit_ch }

// Keep groups with more than 1 replicate, ready for combine replicates
dedup_groupsplit_ch
        .map {
            group, bams -> 
                        if (bams.size() != 1){ //Only keep 'groups' with >1 replicate
                            tuple( groupKey(group, bams.size()), bams)
                        }
        }
        .set{bam_sorted_groups_ch}




// Combine replicates based on 'sample_rep' format, all 'sample' bam files will be merged
process combine_replicates {
    tag "$group"
    publishDir "${params.outdir}/${params.aligner}/samtools/merge/${group}/", pattern:'*', mode: 'copy'

    input:
    tuple val(group), path(bams) from bam_sorted_groups_ch

    output:
    tuple val(group), file("*.merged.bam") into bam_merged_groups_ch

    script:
    """
    samtools merge -n ${group}.merged.bam $bams 
    """
}
  

// Output bam channels into single channel
bam_merged_groups_ch
                .mix(dedup_bam_ch2)
                .into{ dedup_bam_merged_ch; dedup_bam_merged_ch2 }


// samtools index bam files for MethylExtract
process samtools_sort_index {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/samtools/bam/${key}/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from dedup_bam_merged_ch
    
    output:
    tuple val(key), path('*sort.bam') into bam_merged_sorted_ch
    tuple val(key), path('*bai') into bam_merged_indexed_ch

    when:
    params.aligner == 'bwa-meth'

    script:
    """
    samtools sort $bam > ${key}.sort.bam

    samtools index ${key}.sort.bam 
    """
} 

bam_merged_sorted_ch
            .join(bam_merged_indexed_ch)
            .set{bam_sorted_indexed_ch}


// Extract methylation information from Bismark bam files

process methylation_extractor {
    tag "$key"
    publishDir "${params.outdir}/bismark/methylation/$key", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from dedup_bam_merged_ch2

    when:
    params.aligner == 'bismark'

    output:
    path("*txt") into methylation_report_ch
    tuple val(key), path("*.bedGraph.gz") into bedgraph
    tuple val(key), path("*.txt.gz") into methylation_calls
    tuple val(key), path("*.cov.gz") into coverage
    tuple val(key), path("*_splitting_report.txt") into report
    tuple val(key), path("*.M-bias.txt") into mbias
    

    script:
    """
    bismark_methylation_extractor \\
                                    $bam \\
                                    --bedGraph \\
                                    --counts \\
                                    --gzip \\
                                    --report \\
                                    --paired
    """
}


// Extract methylation information from bwa-meth bam files using MethylDackel
// NOTE: requires indexed bam files in same repo as bam files
process methyldackel {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/MethylDackel/$key", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam), path(bai) from bam_sorted_indexed_ch
    path(genome) from params.genome

    output:
    tuple val(key), path("*.bedGraph"), optional: true
    tuple val(key), path("*.methylKit"), optional: true

    when:
    params.aligner == 'bwa-meth'

    script:
    """
    MethylDackel extract \\
                $genome \\
                $bam
    """
}


// Create multiqc report channel based on aligner
if(params.aligner == 'bismark'){

  multiqc_ch = fastqc_ch
                    .mix(fastqc_trimmed_ch)
                    .mix(align_report_ch)                    
                    .mix(dedup_report_ch)
                    .mix(methylation_report_ch)
                    .collect()

}else if(params.aligner == 'bwa-meth'){
  multiqc_ch = fastqc_ch
                    .mix(fastqc_trimmed_ch)
                    .mix(samtools_stats_ch)   
                    .mix(samtools_flag_ch)
                    .mix(dedup_report_ch)
                    .collect()
}

process multiqc {

    publishDir "${params.outdir}/multiqc/", pattern:"multiqc_report.html", mode:'copy'
       
    input:
    path('*') from multiqc_ch
    
    output:
    path('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 




workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong" )
}
