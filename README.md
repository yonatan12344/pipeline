# Reference based viral genome assembly pipeline

# Requirments

1. Nextflow installation

2. Docker Desktop running while running the pipeline

3. Nextflow version 25.10.4.11173 should work

# Pipeline Steps

1. fastq trimming with fastp
2. aggregating post trim statistics with multiqc
3. indexing the user provided reference genome, to generate index files needed for read mapping with bwa-mem2
4. Read mapping with bwa-mem2
5. Read deduplication with samtools
6. Viral Genome assembly with ViralConsesus

* Optional step to view aggregated QC metrics with a subworkflow, to help determine a good phred score cutoff value

# How to Run

1. Open and start Docker desktop in the background

2. Have your raw paired end fastq illumina short reads in a cwd with the main.nf and nextflow.config file.
The reads need to be in a folder called raw_reads, and must match the following regex convention 
'./raw_reads/*_{1,2}.fastq.gz'.

# Output Files
assemblies: contains the final genome assemblies, and tsv & json files detailing the coverage of the assemblies.
deduped_bam: deduplicated bam files
mapped_reads: mapped reads
raw_qc: raw_qc stats
ref_seq_index_files: index files generated for the readmapper to use, based on the ref sequence
trimmed_reads: trimmed fastq files
trimmed_reads_qc: qc information about the post trimmed fastq files

#Quickstart commands

## To view raw read QC info 

nextflow run main.nf -entry assembly --fastq_cutoff 30 --ref_seq RSV_A_ref.fasta 

## To actually perform reference based assembly

nextflow run main.nf -entry assembly --fastq_cutoff 30 --ref_seq RSV_A_ref.fasta 

* replace ref_seq with your chosen reference fasta file
