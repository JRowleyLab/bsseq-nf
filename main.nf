/* 
 * pipeline input parameters 
 */

params.genome = "${baseDir}/../Indexes/bwaIndex/hg38.fa" 
params.samplesheet = "${baseDir}/Samplesheets/samples_test.csv"
params.outdir = "${baseDir}/results"
params.index = "${baseDir}/../Indexes/bwaIndex/"
params.threads = "2"
params.aligner = "bwa-meth"
params.merge = ""  //merge replicates
params.plusmerge = "" //merge replicates but also calculate methylation/duplicates for individual replicates


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
        .into { read_pairs_ch; read_pairs_ch2 }

//  Run fastQC to check quality of reads files

process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", pattern:"{*.html,fastqc_${sample_id}_logs}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from read_pairs_ch

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}


// Download hg38 if no reference is given AND no index is given

process download_ref {

	publishDir "${params.outdir}/reference", pattern: "*.fa", mode:'copy'

	output:
	file("*.fa") into ref_ch

    when: !params.genome & !params.index

	script:
	"""
	wget --no-check-certificate http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_ensembl.dna.fa
	"""
}

genome_ch = params.genome ? Channel.value(file(params.genome)) : ref_ch

// Create index for bismark if none given
// NOT FUNCTIONAL, the output path doesn't pick up the index

process index {

    publishDir "${params.outdir}/${params.aligner}/index", mode: 'copy'
    
    input:
    path genome from genome_ch

    when:
    !params.index
     
    output:
    path 'index' 
    path 'bwameth' into index_out_ch

    script:       
    if( params.aligner == 'bismark' ) 
        """
        bismark_genome_preparation . $genome 
        """
    else if ( params.aligner == 'bwa-meth')
        """
        bwameth.py index $genome
        """
}


index_ch = params.index ? Channel.value(file(params.index)) : index_out_ch



// Trimming reads with Trim Galore

process trimming {

    publishDir "${params.outdir}/trimmed", pattern: "*.fq.gz", mode: 'copy'

    input: 
    tuple val(key), path(reads) from read_pairs_ch2

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

// Merge bam files in a sample 'eg. CONTROL' into a single bam file
// Split bam tuple into groups of sample
if ( params.merge || params.plusmerge) { // 

    // samtools sort bam files, ready for combining replicates (--merge)
    process samtools_sort {
        tag "$key"
        publishDir "${params.outdir}/${params.aligner}/samtools/", pattern:'*', mode: 'copy'

        input:
        tuple val(key), path(bam) from bam_ch
        
        output:
        tuple val(key), path('*sort.bam') into bam_sorted_key_ch

        
        script:
        """
        samtools sort $bam > ${key}.sort.bam
        """
    }


    // TODO: Optimize to move forward with a sample when all replicates have arrived. 
    // Doesn't wait for files to arrive, channel is set as soon as one enters, the channel is then group + all individualss
    // Perhaps using groupKey()
    // bam_sorted_key_ch 
    //                 .map { group_rep, bam ->
    //                                         def(group) = group_rep.split("_")  
    //                                         tuple( group, bam )
    //                                         }
    //                 .groupTuple()
    //                 .set{bam_sorted_groups_ch}


    // Setting groupKey() so that the channel waits for all of a sample to arrive before moving on with the merge
    bam_sorted_key_ch
        .map { group_rep, bam ->
                            def(group) = group_rep.split("_")  
                            tuple( group, bam )
                            }
        .groupTuple()
        .map {
            group, bams -> 
                        tuple( groupKey(group, bams.size()), bams)
        }
        .set{bam_sorted_groups_ch}

//bam_sorted_groups_ch.view()



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
  
}

// Output bam channels into single channel
if(params.merge){
    bam_merged_groups_ch.into{ bam_channel; bam_channel2 }
}else if(params.plusmerge){
    bam_merged_groups_ch.mix(bam_ch2).into{ bam_channel; bam_channel2 }
}else{
    bam_ch2.into{ bam_channel; bam_channel2 }
} 


// samtools index merged (--merge) bam files for MethylExtract
process samtools_sort_index {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/samtools/bam/${key}/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_channel
    
    output:
    tuple val(key), path('*sort.bam') into bam_sorted_channel, bam_sorted_channel2, bam_sorted_channel3
    tuple val(key), path('*bai') into bam_sorted_indexed_channel

    when:
    params.aligner == 'bwa-meth'

    script:
    """
    samtools sort $bam > ${key}.sort.bam

    samtools index ${key}.sort.bam 
    """
} 


// Samtools stats for summary statistics on bwa-meth alignment
process samtools_stat_flagstat {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/samtools/stats/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sorted_channel

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


process samtools_sortn {
        tag "$key"
        publishDir "${params.outdir}/${params.aligner}/samtools/", pattern:'*', mode: 'copy'

        input:
        tuple val(key), path(bam) from bam_channel2
        
        output:
        tuple val(key), path('*sortn.bam') into bam_sortn_channel

        when:
        params.aligner == 'bismark'
        
        script:
        """
        samtools sort -n $bam > ${key}.sortn.bam
        """
    }


// Deduplicate the bam files with bismark
process deduplication_bismark {
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/deduplicate/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sortn_channel

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



// Deduplicate bam files with picard

process deduplication_picard{
    tag "$key"
    publishDir "${params.outdir}/${params.aligner}/picard/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sorted_channel3 
    path(genome) from params.genome

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
    dedup_bam_bismark_ch.into{ dedup_bam_ch; dedup_bam_ch2 }
    dedup_report_ch = dedup_report_bismark_ch
}else{
    dedup_bam_picard_ch.into{ dedup_bam_ch; dedup_bam_ch2 }
    dedup_report_ch = dedup_report_picard_ch
} 
    

// Extract methylation information from Bismark bam files

process methylation_extractor {
    tag "$key"
    publishDir "${params.outdir}/bismark/methylation/$key", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from dedup_bam_ch

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
    tuple val(key), path(bam) from dedup_bam_ch2
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
