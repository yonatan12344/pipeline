params.fastq_cutoff = 20
params.ref_seq = 'please provide ref sequence'

process fastqc_view {
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "results/raw_qc", mode: 'copy'
    
    input:
    path reads

    output:
    path '*.html', emit: fastqc_html
    path '*.zip', emit: fastqc_zip

    script:
    """
    fastqc $reads
    """
}

process multiqc_raw {
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "results/raw_qc", mode: 'copy'

    input:
    path reads

    output:
    path 'multiqc_report.html', emit: fastqc_html

    script:
    """
    multiqc $reads
    """
}

process fastp_trimming {
    container 'quay.io/biocontainers/fastp:1.3.3--h43da1c4_0'
    publishDir "results/trimmed_reads", mode: 'copy'

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("cleaned*.fastq.gz"), emit: reads
    path '*.json', emit: jsons

    script:
    """
    fastp -i ${reads[0]}  \\
    -I ${reads[1]}  \\
    -o cleaned_${reads[0].simpleName}.fastq.gz \\
    -O cleaned_${reads[1].simpleName}.fastq.gz \\
    -q ${params.fastq_cutoff} \\
    --detect_adapter_for_pe \\
    -h ${sampleId}.html \\
    -j ${sampleId}.json
    """
}

process post_trim_multiqc {
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "results/trimmed_reads_qc", mode: 'copy'

    input:
    path json_files

    output:
    path "multiqc_report.html", emit: fastp_html
    path "multiqc_data", emit: fastp_data

    script:
    """
    multiqc .
    """
}

process ref_index_files {
    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0'
    publishDir "results/ref_seq_index_files"

    input:
    path ref_seq

    output:
    path "${ref_seq}.*", emit: index_files

    script:
    """
    bwa-mem2 index "${ref_seq}"
    """
}

process read_mapping {
    container 'quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0'
    publishDir "results/mapped_reads"

    input:
    tuple val(sampleId), path(reads)
    path ref
    path index_files

    output:
    tuple val(sampleId), path("${sampleId}.sam"), emit: sam_file

    script:
    """
    bwa-mem2 mem -t 2 $ref ${reads[0]} ${reads[1]} > ${sampleId}.sam
    """
}

process dedup_bam {
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'
    publishDir "results/deduped_bam", mode: 'copy'

    input:
    tuple val(sampleId), path(sam_file)

    output:
    tuple val(sampleId), path("${sampleId}.dedup.bam"), emit: bam_file

    script:
    """
    # 1. Convert to BAM and sort by name (required for fixmate)
    samtools sort -n -o name_sorted.bam ${sam_file}

    # 2. Add MS tags required by markdup
    samtools fixmate -m name_sorted.bam fixmate.bam

    # 3. Sort by coordinate (required for markdup)
    samtools sort -o coord_sorted.bam fixmate.bam

    # 4. Remove duplicates
    samtools markdup -r coord_sorted.bam ${sampleId}.dedup.bam

    # Clean up intermediate files
    rm name_sorted.bam fixmate.bam coord_sorted.bam
    """
}

process viral_consesus {
    container 'niemasd/viral_consensus:0.0.5'
    publishDir "results/assemblies"

    input:
    tuple val(sampleId), path(sam_file)
    path ref

    output:
    path("*.fasta"), emit: fasta_assembly
    path("*.tsv"), emit: tsv_file
    path("*.json"), emit: json_file 

    script:
    """
    viral_consensus -i $sam_file -r $ref -o ${sampleId}.fasta -op ${sampleId}.tsv -oi ${sampleId}.json
    """
}

workflow raw_read_view {
    raw_reads = channel.fromPath("./raw_reads/*.fastq.gz")
    fastqc_view(raw_reads)
    multiqc_raw(fastqc_view.out.fastqc_zip.collect())
}

workflow assembly {
    raw_reads = Channel.fromFilePairs('./raw_reads/*_{1,2}.fastq.gz')

    fastp_out = fastp_trimming(raw_reads)
    
    post_trim_multiqc(fastp_out.jsons.collect())

    ref = Channel.value(file(params.ref_seq))

    ref_index = ref_index_files(ref)

    read_mapping(fastp_out.reads, ref, ref_index.index_files.collect())
    
    dedup_bam(read_mapping.out.sam_file)

    viral_consesus(dedup_bam.out.bam_file, ref)
}
