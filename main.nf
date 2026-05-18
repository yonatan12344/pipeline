nextflow.enable.dsl = 2

// fastqc and multiqc to view initial quality
process raw_qc_view {
	conda 'fastqc multiqc'
	publishDir 'raw_qc_metrics', mode: 'copy'
input:
	path "./raw_reads/*.fastq.gz" 

output:
	path "*.html", emit: html_files
	path "*.zip", emit: zip_files

script:

"""
mkdir -p raw_reads_qc
fastqc ./raw_reads/*.fastq* -o ./raw_reads_qc
multiqc ./raw_reads -o ./raw_reads_qc
"""

}

workflow {
raw_fasta = channel.fromFilePairs('./raw_reads/*{1,2}.fastq.gz')
raw_qc_view(raw_fasta)
    }