process extract_rnd_barcode {
	tag {query.simpleName}

	input:
	path query // from fastq_quality_filter_to_barcode_extraction

	output:
	path "${query.baseName}.rndBarcode.fastq" , emit: rnd_barcode// into fastq_barcode_extraction_to_split_exp_barcode
	path ("${query.simpleName}.log", emit: barcode_log) // into collect_statistics_rnd_barcode_extraction

	"""
	umi_tools extract --stdin ${query} --bc-pattern ${params.barcode_pattern} --log ${query.simpleName}.log --stdout ${query.baseName}.rndBarcode.fastq
	"""
}

process check_barcode_file {

	input:
	path "barcodes" // from barcode_file

	output:
	path "checkedBarcodes" // into checked_barcodes

	"""
	check_barcode_file.py barcodes > checkedBarcodes
	"""
}

//FASTX Toolkit 0.0.14
process split_exp_barcode {
	tag {query.simpleName}

	input:
	tuple file(query), file(barcodes) // from fastq_barcode_extraction_to_split_exp_barcode.combine(checked_barcodes)

	output:
	path "barcode_split.log" , emit: log  // into (collect_statistics_split,log_split_experimental_barcode)
	path "output/*.fastq" , emit: fastq_split_barcode_to_remove_exp_barcode // into fastq_split_barcode_to_remove_exp_barcode

	script:
	if(params.barcode_mismatches == 0)
		"""
		mkdir output
		cat ${query} | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --exact --prefix ./output/ --suffix .fastq > log_split.txt
		"""
	else
		"""
		mkdir output
		cat ${query} | fastx_barcode_splitter.pl --bcfile ${barcodes} --bol --mismatches ${params.barcode_mismatches} --prefix ./output/ --suffix .fastq > barcode_split.log
		"""
}

process get_length_exp_barcode {

	input:
	val(pattern) // from val_barcode_pattern

	output:
	env(exp_barcode_length) // into int_length_exp_barcode

	"""
	exp_barcode_length=\$((\$(echo -n ${pattern} | sed 's/[^Xx]//g' | wc -m) + 1)) 
	"""
}

//FASTX Toolkit 0.0.14
process remove_exp_barcode {
	tag {query.simpleName}

	input:
	tuple file(query), val(length_exp_barcode) // from fastq_split_barcode_to_remove_exp_barcode.flatten().filter{ it.size() > 0 }.combine(int_length_exp_barcode)

	output:
	path "*.preprocessed.fastq" // into fastq_remove_barcode_to_collect

	"""
	fastx_trimmer -f ${length_exp_barcode} -i ${query} -o ${query.baseName}.preprocessed.fastq
	"""
}



process merge_preprocessed_reads {
	tag {name}

	input:
	tuple val(name), file("query") // from fastq_collect_preprocessed_to_merge

	output:
	file("${name}.fastq") // into fastq_merge_preprocessed_to_alignment

	"""
	cat ${query} > ${name}.fastq
	"""
}

process generate_barcode_barplot {
	publishDir "${params.output}/statistics", mode: 'copy'

	input:
	path(query) // from log_split_experimental_barcode.first()

	output:
	path("${query.baseName}.png") // into output_experimental_barcode_distribution_png

	"""
	plot_experimental_barcode_distribution.R --logs ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
	"""
}