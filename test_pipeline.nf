
nextflow.preview.dsl=2


//FastQC v0.11.9
process quality_control {
	tag {query.simpleName}
	
	input:
	file query // from input_reads_QC

	output:
	file "${query.baseName}*" // intoout1, collect_statistics_qc1

	"""
	fastqc ${query} -o .
	"""
}


