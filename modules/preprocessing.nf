// this files conains the processes for the preprocessing step of the PARANOID pipeline

//FastQC v0.11.9
process quality_control {
	tag {query.simpleName}
	
	input:
	path query // from input_reads_QC

	output:
	path "${query.baseName}*" // into out1, collect_statistics_qc1

	"""
	fastqc ${query} -o .
	"""
}

//FastQC v0.11.9
process quality_control_2 {
	tag {query.simpleName}
	
	input:
	path query // from fastq_quality_filter_to_quality_control_2.first()

	output:
	path "quality-control-2*" // into qc_2_out, collect_statistics_qc2

	"""
	cat ${query} > quality-control-2.fastq
	fastqc quality-control-2.fastq -o .
	"""
}


process adapter_removal {
	tag {query.simpleName}

	input:
	path query // from input_reads_processing

	output:
	path "${query}_trimmed.fq", emit: trimmed_fq // into reads_qualityFilter  
	path "${query}_trimming_report.txt", emit: trimming_report // into collect_statistics_adapter

	"""
	trim_galore --cores ${task.cpus} --basename ${query} -o . --length ${params.min_length} ${query} --quality 0
	"""
}

//FASTX Toolkit 0.0.14
process quality_filter {
	tag {query.simpleName}

	input:
	path query // from reads_qualityFilter

	output:
	path "${query.baseName}.qual-filter.fastq" , emit: filter_fq// into fastq_quality_filter_to_barcode_extraction, fastq_quality_filter_to_quality_control_2
	path 'summary-quality-filter.txt' , emit: filter_report // into collect_statistics_quality_filter


	"""
	fastq_quality_filter -v -q ${params.min_qual} -p ${params.min_percent_qual_filter} -i ${query} -o ${query.baseName}.qual-filter.fastq > summary-quality-filter.txt
	"""
}


