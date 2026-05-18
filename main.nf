params.ref_seq = null
process read_view {
	conda 'fastp multiqc'
	publishDir Results, mode: 'copy'
	script:
	"""
	"""
}


process read_mapping {
	conda 'fastp multiqc'
	publishDir Results, mode: 'copy'
	script:
	"""
	"""
}

process assembly {
	
}

 