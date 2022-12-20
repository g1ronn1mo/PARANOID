#!/usr/bin/env nextflow


//FastQC v0.11.9
process quality_control {
	tag {query.simpleName}
	
	input:
	file query from input_reads_QC

	output:
	file "${query.baseName}*" into out1, collect_statistics_qc1

	"""
	fastqc ${query} -o .
	"""
}

process adapter_removal {
	tag {query.simpleName}

	input:
	file query from input_reads_processing

	output:
	file "${query}_trimmed.fq" into reads_qualityFilter
	file "${query}_trimming_report.txt" into collect_statistics_adapter

	"""
	trim_galore --cores ${task.cpus} --basename ${query} -o . --length ${params.min_length} ${query} --quality 0
	"""
}

//FASTX Toolkit 0.0.14
process quality_filter {
	tag {query.simpleName}

	input:
	file query from reads_qualityFilter

	output:
	file "${query.baseName}.qual-filter.fastq" into fastq_quality_filter_to_barcode_extraction, fastq_quality_filter_to_quality_control_2
	file 'summary-quality-filter.txt' into collect_statistics_quality_filter


	"""
	fastq_quality_filter -v -q ${params.min_qual} -p ${params.min_percent_qual_filter} -i ${query} -o ${query.baseName}.qual-filter.fastq > summary-quality-filter.txt
	"""
}

//FastQC v0.11.9
process quality_control_2 {
	tag {query.simpleName}
	
	input:
	file query from fastq_quality_filter_to_quality_control_2.first()

	output:
	file "quality-control-2*" into qc_2_out, collect_statistics_qc2

	"""
	cat ${query} > quality-control-2.fastq
	fastqc quality-control-2.fastq -o .
	"""
}

process extract_rnd_barcode {
	tag {query.simpleName}

	input:
	file query from fastq_quality_filter_to_barcode_extraction

	output:
	file "${query.baseName}.rndBarcode.fastq" into fastq_barcode_extraction_to_split_exp_barcode
	file("${query.simpleName}.log") into collect_statistics_rnd_barcode_extraction

	"""
	umi_tools extract --stdin ${query} --bc-pattern ${params.barcode_pattern} --log ${query.simpleName}.log --stdout ${query.baseName}.rndBarcode.fastq
	"""
}

process check_barcode_file {

	input:
	file "barcodes" from barcode_file

	output:
	file "checkedBarcodes" into checked_barcodes

	"""
	check_barcode_file.py barcodes > checkedBarcodes
	"""
}

//FASTX Toolkit 0.0.14
process split_exp_barcode {
	tag {query.simpleName}

	input:
	set file(query), file(barcodes) from fastq_barcode_extraction_to_split_exp_barcode.combine(checked_barcodes)

	output:
	file "barcode_split.log" into (collect_statistics_split,log_split_experimental_barcode)
	file "output/*.fastq" into fastq_split_barcode_to_remove_exp_barcode

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
	val(pattern) from val_barcode_pattern

	output:
	env(exp_barcode_length) into int_length_exp_barcode

	"""
	exp_barcode_length=\$((\$(echo -n ${pattern} | sed 's/[^Xx]//g' | wc -m) + 1)) 
	"""
}

//FASTX Toolkit 0.0.14
process remove_exp_barcode {
	tag {query.simpleName}

	input:
	set file(query), val(length_exp_barcode) from fastq_split_barcode_to_remove_exp_barcode.flatten().filter{ it.size() > 0 }.combine(int_length_exp_barcode)

	output:
	file "*.preprocessed.fastq" into fastq_remove_barcode_to_collect

	"""
	fastx_trimmer -f ${length_exp_barcode} -i ${query} -o ${query.baseName}.preprocessed.fastq
	"""
}

fastq_remove_barcode_to_collect
	.map{file -> tuple(file.name - ~/\.[\w.]+.fastq$/,file)}
	.groupTuple()
	.set{ fastq_collect_preprocessed_to_merge }

process merge_preprocessed_reads {
	tag {name}

	input:
	set val(name), file("query") from fastq_collect_preprocessed_to_merge

	output:
	file("${name}.fastq") into fastq_merge_preprocessed_to_alignment

	"""
	cat ${query} > ${name}.fastq
	"""
}

if ( params.domain == 'pro' || params.map_to_transcripts == true){
	//bowtie2 version 2.3.5.1
	process build_index_bowtie {

		input:
		file ref from reference_to_mapping

		output:
		set file("${ref}"), file("${ref}.*") into bowtie_index_build_to_mapping

		"""
		bowtie2-build ${ref} ${ref}
		"""
	}

	process mapping_bowtie{
		tag {query.simpleName}

		input:
		set file(ref), file(index) from bowtie_index_build_to_mapping.first()
		file query from fastq_merge_preprocessed_to_alignment

		output:
		file "${query.baseName}.bam" into (bam_mapping_to_filter_empty, bam_mapping_to_sort)
		file "${query.simpleName}.statistics.txt" into collect_statistics_mapping

		"""
		bowtie2 --no-unal -q -p ${task.cpus} -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
		"""
	}
} else if ( params.domain == 'eu' ) {
	//Version 2.7.3a
	process build_index_STAR {

		input:
		file referenceGenome from reference_to_mapping
		file gtf from sjdbGTFfile

		output:
		file index into star_index_build_to_mapping

		script:
		if(params.annotation == 'NO_FILE')
			"""
			mkdir index
			STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} 
			"""
		else
			"""
			mkdir index
			STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} --sjdbGTFfile ${gtf}
			"""
	}

	process mapping_STAR{
		tag {query.simpleName}

		input:
		set file(query), file(indexDir) from fastq_merge_preprocessed_to_alignment.combine(star_index_build_to_mapping)

		output:
		file("${query.baseName}.Aligned.sortedByCoord.out.bam") into bam_mapping_to_filter_empty
		file("${query.baseName}.Log.*") into collect_statistics_mapping

		"""
		STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --alignEndsType Extend5pOfRead1 --outSAMtype BAM SortedByCoordinate
		"""
	}
}

process filter_empty_bams{
	tag {query.simpleName}

	input:
	file query from bam_mapping_to_filter_empty

	output:
	file "${query.baseName}.filtered.bam" optional true into bam_filter_empty_to_split
	file "${query.simpleName}.no_alignments.txt" optional true into log_experiments_without_alignments

	"""
	if [[ \$(samtools view ${query} | wc -l) == 0 ]]; then
		echo "${query.simpleName} contains no alignable reads" > ${query.simpleName}.no_alignments.txt
	else
		mv ${query} ${query.baseName}.filtered.bam
	fi
	"""
}

if( params.map_to_transcripts == false && params.speed ) {
	process split_bam_by_chromosome{
		tag {query.simpleName}

		input:
		file(query) from bam_filter_empty_to_split

		output:
		file("*.bam") into bam_split_to_sort

		"""
		bamtools split -in ${query} -reference
		"""
	}
} else {
	bam_filter_empty_to_split.set{ bam_split_to_sort }
}

process sort_bam{
	tag {query.simpleName}

	input:
	file query from bam_split_to_sort.flatten()

	output:
	set file("${query.baseName}.sorted.bam"), file("${query.baseName}.sorted.bam.bai")  into bam_sort_to_deduplicate

	"""
	samtools sort ${query} -o ${query.baseName}.sorted.bam
	samtools index ${query.baseName}.sorted.bam
	"""
}

process deduplicate{
	publishDir "${params.output}/statistics/PCR-deduplication", mode: 'copy', pattern: "${query.baseName}.deduplicated.log*"
	tag {query.simpleName}
	memory { 20.GB + 1.B * query.size() }

	input:
	set file(query), file(index) from bam_sort_to_deduplicate

	output:
	file "${query.baseName}.deduplicated.bam" into (bam_deduplicate_to_sort,bam_deduplicate_to_index)
	file "${query.baseName}.deduplicated.log*" into (log_deduplicate_to_collect_statistics,log_deduplicate_to_output)

	"""
	umi_tools dedup --random-seed=42 -I ${query} --output-stats ${query.baseName}.deduplicated.log -S ${query.baseName}.deduplicated.bam
	"""
}

bam_deduplicate_to_sort
	.map{file -> tuple(file.name - ~/\.[\w.]+.bam$/,file)}
	.groupTuple()
	.set{ bam_dedup_sort_to_merge }

process merge_deduplicated_bam {
	tag {name}

	input:
	set val(name), file(query) from bam_dedup_sort_to_merge

	output:
	file("${name}.bam") into (bam_merge_to_calculate_crosslinks, bam_merge_to_extract_transcripts, bam_merge_to_pureCLIP)
	file("${name}.bam") into bam_alignment_to_sort

	"""
	samtools merge -c -p ${name}.bam ${query}
	"""
}

process sort_and_index_alignment{
	tag {query.simpleName}
	publishDir "${params.output}/alignments", mode: 'copy', pattern: "${query.simpleName}.sorted.bam*"

	input:
	file(query) from bam_alignment_to_sort

	output:
	file("${query.simpleName}.sorted.bam") into bam_sorted_to_igv_session
	file("${query.simpleName}.sorted.bam*")

	"""
	samtools sort ${query} > ${query.simpleName}.sorted.bam
	samtools index ${query.simpleName}.sorted.bam
	"""
}

if (params.map_to_transcripts == true){
	process count_hits {
		tag {bam.simpleName}

		publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

		input:
		file bam from bam_merge_to_extract_transcripts

		output:
		file("${bam.baseName}.hits.tsv") into tsv_count_to_get_hits

		"""
		samtools view ${bam} | cut -f3 | sort | uniq -c | sort -nr > ${bam.baseName}.hits.tsv
		"""
	}

	process get_top_hits {

		publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

		input:
		file tsv from tsv_count_to_get_hits.flatten().toList()

		output:
		file("transcript-targets-top${params.number_top_transcripts}.txt") into (txt_get_hits_to_extract_alignments,txt_get_hits_to_extract_sequences)

		"""
		head -${params.number_top_transcripts} -q *.tsv  | rev | cut -f1 -d' ' | rev | sort | uniq > transcript-targets-top${params.number_top_transcripts}.txt
		"""
	}

	process index_alignments {
		tag {bam.simpleName}

		input:
		file bam from bam_deduplicate_to_index

		output:
		set file("${bam.baseName}.sorted.bam"), file("${bam.baseName}.sorted.bam.bai") into bam_index_to_extract_alignments

		"""
		samtools sort ${bam} -o ${bam.baseName}.sorted.bam
		samtools index ${bam.baseName}.sorted.bam
		"""
	}

	process extract_top_alignments {
		tag {bam.simpleName}

		input:
		set file(bam), file(bai), file(txt_alignments) from bam_index_to_extract_alignments.combine(txt_get_hits_to_extract_alignments)

		output:
		file("${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam") into bam_extract_alignments_to_calc_crosslink

		"""
		samtools view -hb ${bam} `cat ${txt_alignments}` > ${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam
		"""
	}

	process remove_newlines {

		input:
		file ref from reference_to_extract_transcripts

		output:
		file("${ref.baseName}.removed_newlines.fna") into fasta_rm_newline_to_extract_sequences

		"""
		awk '!/^>/ {printf "%s", \$0; n = "\\n" } /^>/ { print n \$0; n = "" } END { printf "%s", n }' ${ref} > ${ref.baseName}.removed_newlines.fna
		"""
	}

	process extract_top_transcript_sequences {
		tag {txt_sequences.simpleName}

		publishDir "${params.output}/transcripts", mode: 'copy'

		input:
		set file(txt_sequences), file(ref) from txt_get_hits_to_extract_sequences.combine(fasta_rm_newline_to_extract_sequences)

		output:
		file("${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna") into (top_transcripts_to_collect,top_transcripts_to_strand_preference)


		"""
		egrep -A1 --no-group-separator -f ${txt_sequences} ${ref} > ${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna
		"""
	}

	bam_extract_alignments_to_calc_crosslink
	.into{collected_bam_files; collected_bam_files_to_sort}
} else {
	bam_merge_to_calculate_crosslinks
	.into{collected_bam_files; collected_bam_files_to_sort}
}

process sort_bam_before_strand_pref {
	tag {query.baseName}

	input:
	file(query) from collected_bam_files_to_sort

	output:
	file("${query.baseName}.sorted.bam") into bam_sort_to_group

	"""
	samtools sort ${query} > ${query.baseName}.sorted.bam
	"""
}

if(params.merge_replicates == true){
	bam_sort_to_group
	.map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*(.sorted)?.bam$/,file)} 
	.groupTuple()
	.set{grouped_bam_to_strand_preference}
} else {
	bam_sort_to_group
	.map{file -> tuple(file.name - ~/(_filtered_top\d*)?(.sorted)?.bam$/,file)}
	.set{grouped_bam_to_strand_preference}
}

if (params.map_to_transcripts == true){
	top_transcripts_to_strand_preference
	.set{choose_correct_reference_file}
} else {
	reference_to_strand_preference
	.set{choose_correct_reference_file}
}

process determine_strand_preference {
	tag {name}
	publishDir "${params.output}/strand-distribution", mode: 'copy', pattern: "${name}.strand_proportion.txt"

	input:
	set val(name),file(query),file(reference) from grouped_bam_to_strand_preference.combine(choose_correct_reference_file)

	output:
	file("${name}.strand_proportion.txt") into txt_determine_strand_to_visualize

	"""
	egrep '^>' ${reference} | cut -f1 -d' ' | cut -c2- > references.txt
	for i in ${query}
	do
		samtools index \$i
	done
	touch ${name}.strand_proportion.txt
	echo -e "chromosome\tforward\treverse" >> ${name}.strand_proportion.txt
	while read r; do
  		echo -e "\$r\t\$(for i in *.sorted.bam; do samtools view -F 20 -q ${params.mapq} \$i \$r | wc -l; done | paste -s -d+ | bc)\t\$(for i in *.sorted.bam; do samtools view -f 16 -q ${params.mapq} \$i \$r | wc -l; done | paste -s -d+ | bc)" >> ${name}.strand_proportion.txt;
	done <references.txt
	"""
}

process visualize_strand_preference {
	publishDir "${params.output}/strand-distribution/visualization", mode: 'copy', pattern: "${strand.simpleName}.png"

	input:
	file(strand) from txt_determine_strand_to_visualize

	output:
	file("${strand.simpleName}.png")

	"""
	visualize_strand_distribution.R --input ${strand} --output ${strand.simpleName} --type png
	"""
}

process get_chromosome_sizes{
	input:
	file(ref) from reference_to_chrom_sizes

	output:
	file("${ref.simpleName}.chromosome_sizes.txt") into (chrom_sizes_to_cross_link_calculation,chrom_sizes_to_bigWig)

	"""
	samtools faidx ${ref}
	cut -f1,2 ${ref}.fai > ${ref.simpleName}.chromosome_sizes.txt
	"""
}

process calculate_crosslink_sites{
	tag {query.simpleName}
	publishDir "${params.output}/cross-link-sites/wig", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

	input:
	set file(query), file(chrom_sizes) from collected_bam_files.combine(chrom_sizes_to_cross_link_calculation)

	output:
	file "${query.simpleName}.wig2" optional true into (wig_calculate_crosslink_to_group_samples,wig2_calc_cl_sites_to_split)
	set val("cross-link-sites"), file("${query.simpleName}_forward.wig"), file("${query.simpleName}_reverse.wig") optional true into wig_cross_link_sites_to_transform
	file "${query.simpleName}_{forward,reverse}.wig" optional true into wig_calculate_cl_sites_to_split_strand

	"""
	create-wig-from-bam.py --input ${query} --mapq ${params.mapq} --chrom_sizes ${chrom_sizes} --output ${query.simpleName}.wig2
	if [[ -f "${query.simpleName}.wig2" ]]; then
		wig2-to-wig.py --input ${query.simpleName}.wig2 --output ${query.simpleName}
	fi
	"""
}

//From here on only further analyses

if( params.merge_replicates == true ){

	//groups files according to their experiment
	wig_calculate_crosslink_to_group_samples
	.map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*.wig2$/,file)} 
	.groupTuple()
	.set{wig2_grouped_samples_to_merge}

	process split_wig2_for_correlation{
		tag {query.simpleName}

		input:
		file(query) from wig2_calc_cl_sites_to_split

		output:
		file("${query.simpleName}_forward.wig") optional true into split_wig_forward_to_correlation
		file("${query.simpleName}_reverse.wig") optional true into split_wig_reverse_to_correlation

		"""
		wig2-to-wig.py --input ${query} --output ${query.simpleName}
		"""
	}

	split_wig_forward_to_correlation
	.map{file -> tuple(file.name - ~/_rep_\d*_forward.wig$/,file)} 
	.groupTuple()
	.set{group_forward}

	/*process calc_wig_correlation{
		tag {name}
		echo true
		cache false

		input:
		set val(name),file(query) from group_forward

		output:
		val("${name}")

		when:
		{query}.size() >= 2

		"""
		echo ${query}
		"""
	}*/

	process merge_wigs{
		tag {name}
		publishDir "${params.output}/cross-link-sites-merged/wig", mode: 'copy', pattern: "${name}_{forward,reverse}.wig"

		input:
		set name, file(query) from wig2_grouped_samples_to_merge

		output:
		file "${name}.wig2" into collected_wig_files
		set val("cross-link-sites-merged"), file("${name}_forward.wig"), file("${name}_reverse.wig") into wig_merged_cross_link_sites_to_transform
		file "${name}_{forward,reverse}.wig" into output

		script:
		if(name !=~ "unmatched*")
			"""
			merge-wig.py --wig ${query} --output ${name}.wig2
			wig2-to-wig.py --input ${name}.wig2 --output ${name}
			"""
		else
			"""
			merge-wig.py --wig ${query} --output unmatched.wig2
			wig2-to-wig.py --input unmatched.wig2 --output unmatched
			"""
	}
} else {
	wig_calculate_crosslink_to_group_samples
		.set{collected_wig_files}
}

// Transformation of cross-link sites to different formats
if( params.merge_replicates == true ) {
	wig_cross_link_sites_to_transform.mix(wig_merged_cross_link_sites_to_transform)
	.set{wig_cross_link_sites_to_bigWig}
} else {
	wig_cross_link_sites_to_transform
	.set{wig_cross_link_sites_to_bigWig}
}

process wig_to_bigWig{
	tag {forward.simpleName}
	publishDir "${params.output}/${out_dir}/bigWig", mode: 'copy', pattern: "*.bw"

	input:
	set val(out_dir), file(forward), file(reverse), file(chrom_sizes) from wig_cross_link_sites_to_bigWig.combine(chrom_sizes_to_bigWig)

	output:
	file("*.bw") optional true
	set val(out_dir), file("*.bw") optional true into big_wig_to_igv_session
	set val(out_dir), file("${reverse.baseName}.bw") optional true into big_wig_reverse_to_convert_to_bedgraph
	set val(out_dir), file("${forward.baseName}.bw") optional true into big_wig_forward_to_convert_to_bedgraph

	"""
	if [[ \$(cat ${forward} | wc -l) > 1 ]]; then
		wigToBigWig ${forward} ${chrom_sizes} ${forward.baseName}.bw
	fi
	if [[ \$(cat ${reverse} | wc -l) > 1 ]]; then
		wigToBigWig ${reverse} ${chrom_sizes} ${reverse.baseName}.bw
	fi
	"""
}

process bigWig_to_bedgraph{
	tag {bigWig.simpleName}
	publishDir "${params.output}/${out_dir}/bedgraph", mode: 'copy', pattern: "*.bedgraph"

	input:
	set val(out_dir), file(bigWig) from big_wig_forward_to_convert_to_bedgraph.mix(big_wig_reverse_to_convert_to_bedgraph)

	output:
	file("*.bedgraph")

	"""
	bigWigToBedGraph ${bigWig} ${bigWig.baseName}.bedgraph
	"""
}

// Generate one channel per postprocessing analysis
collected_wig_files.into{ collected_wig_2_to_RNA_subtypes_distribution; collected_wig_2_to_sequence_extraction; collected_wig_2_to_peak_distance; collect_wig_2_to_peak_height_histogram }

if (params.omit_peak_calling == false){
	process index_for_peak_calling {
		tag{query.simpleName}

		input:
		file(query) from bam_merge_to_pureCLIP

		output:
		set file("${query.simpleName}.sorted.bam"), file("${query.simpleName}.sorted.bam.bai") into bambai_index_to_peak_calling

		"""
		samtools sort ${query} > ${query.simpleName}.sorted.bam
		samtools index ${query.simpleName}.sorted.bam
		"""
	}

	process pureCLIP {
		tag{bam.simpleName}
		errorStrategy 'ignore' //TODO: is supposed to be only temporal. Need to find a solution for: ERROR: Emission probability became 0.0! This might be due to artifacts or outliers.

		publishDir "${params.output}/peak_calling", mode: 'copy', pattern: "${bam.simpleName}.pureCLIP_crosslink_{sites,regions}.bed"

		input:
		set file(bam), file(bai), file(ref) from bambai_index_to_peak_calling.combine(reference_to_pureCLIP)

		output:
		file("${bam.simpleName}.pureCLIP_crosslink_{sites,regions}.bed")
		file("${bam.simpleName}.pureCLIP_crosslink_sites.bed") into bed_peak_calling_to_sequence_extraction
		file("${bam.simpleName}.pureCLIP_crosslink_sites.params") into params_peak_calling_to_collect_statistics

		script:
		if(params.peak_calling_for_high_coverage == true && params.peak_calling_regions == true)
			"""
			pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -or ${bam.simpleName}.pureCLIP_crosslink_regions.bed -dm ${params.peak_calling_regions_width} -mtc 5000 -mtc2 5000 -ld
			"""
		else if(params.peak_calling_for_high_coverage == true && params.peak_calling_regions == false)
			"""
			pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -mtc 5000 -mtc2 5000 -ld
			"""
		else if(params.peak_calling_for_high_coverage == false && params.peak_calling_regions == true)
			"""
			pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed -or ${bam.simpleName}.pureCLIP_crosslink_regions.bed -dm ${params.peak_calling_regions_width}
			"""
		else if(params.peak_calling_for_high_coverage == false && params.peak_calling_regions == false)
			"""
			pureclip -i ${bam} -bai ${bai} -g ${ref} -nt ${task.cpus} -o ${bam.simpleName}.pureCLIP_crosslink_sites.bed
			"""
	}
}

process split_wig_2_for_peak_height_hist {
	tag {query.simpleName}

	input:
	file(query) from collect_wig_2_to_peak_height_histogram

	output:
	set val("${query.simpleName}"), file("${query.simpleName}_forward.wig"), file("${query.simpleName}_reverse.wig") into split_wig2_to_generate_peak_height_histogram

	"""
	wig2-to-wig.py --input ${query} --output ${query.simpleName}
	"""
}

process generate_peak_height_histogram {
	tag {query}
	publishDir "${params.output}/peak_height_distribution", mode: 'copy'

	input:
	set val(query), file(forward), file(reverse) from split_wig2_to_generate_peak_height_histogram

	output:
	file("${query}.png")

	"""
	generate-peak-height-histogram.R --input . --output ${query} --type png --color "${params.color_barplot}" --percentile ${params.percentile}
	"""
}

if ( params.annotation != 'NO_FILE'){
	process wig_to_bam {
		tag {query.simpleName}

		input:
		file(query) from collected_wig_2_to_RNA_subtypes_distribution

		output:
		file("${query.baseName}.bam") into bam_convert_to_feature_counts

		"""
		wig-to-bam.py --input ${query} --output ${query.baseName}.bam
		"""
	}

	process feature_counts {
		tag {query.simpleName}

		input:
		set file(query), val(rna_subtypes), file(annotation) from bam_convert_to_feature_counts.combine(rna_subtypes_to_feature_counts).combine(annotation_to_RNA_subtypes_distribution)

		output:
		file("${query.simpleName}.${rna_subtypes}.tsv") into tsv_feature_counts_to_sort


		"""
		featureCounts -T ${task.cpus} -t ${rna_subtypes} -g ${params.gene_id} -a ${annotation} -R CORE -M -o ${query.simpleName} ${query}
		mv ${query}.featureCounts ${query.simpleName}.${rna_subtypes}.tsv
		"""
	}

	tsv_feature_counts_to_sort
		.map{file -> tuple(file.name - ~/\.[\w.]+.tsv$/,file)}
		.groupTuple()
		.set{tsv_sort_to_calculate_distribution}

	process get_RNA_subtypes_distribution {
		tag {name}

		publishDir "${params.output}/RNA_subtypes", mode: 'copy'

		input:
		val(subtypes) from rna_subtypes_to_distribution.collect()
		set val(name), file(query) from tsv_sort_to_calculate_distribution

		output:
		file("${name}.subtypes_distribution.tsv") into (output_rna_subtypes_tsv,tsv_rna_subtypes_distribution_to_barplot)
		file("${name}.subtype.log") optional true into log_subtype_warnings

		script:
		subtypes_as_string = subtypes.join(' ')
		"""
		calc-RNA-subtypes-distribution.py --input ${query} --rna_subtypes ${subtypes_as_string} --output ${name}.subtypes_distribution.tsv > ${name}.subtype.log
		"""
	}

	process generate_RNA_subtypes_barplot {
		publishDir "${params.output}/RNA_subtypes", mode: 'copy'

		input:
		file(query) from tsv_rna_subtypes_distribution_to_barplot

		output:
		file("${query.baseName}.png") into output_rna_subtypes_png

		"""
		RNA_subtypes_barcharts.R --input ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
		"""
	}

	process collect_subtype_analysis_errors {
		publishDir "${params.output}/statistics", mode: 'copy', pattern: 'subtype-analysis-warnings.txt'

		input:
		file(query) from log_subtype_warnings.flatten().toList()

		output:
		file("subtype-analysis-warnings.txt") optional true into output_log_subtype_analysis_warnings

		"""
		for i in ${query}
		do
			if [[ ! \$(cat \$i | wc -l) == 0 ]]; then
				echo "File: \$i" >> subtype-analysis-warnings.txt
				cat \$i >> subtype-analysis-warnings.txt
			fi
		done
		"""
	}
}

if (params.omit_sequence_extraction == false) {

	if(params.omit_peak_calling == false){
		bed_peak_calling_to_sequence_extraction.set{peaks_for_sequence_extraction}
	} else {
		collected_wig_2_to_sequence_extraction.set{peaks_for_sequence_extraction}
	}

	process sequence_extraction {
		tag {query.simpleName}

		publishDir "${params.output}/extracted_sequences", mode: 'copy', pattern: "*.extracted-sequences.*"

		input:
		set file(query),file(reference) from peaks_for_sequence_extraction.combine(reference_to_extract_sequences)

		output:
		file "*.extracted-sequences.fasta" optional true into extracted_sequences
		file "*.extracted-sequences.*" optional true into extracted_sequences_output
		
		script:
		if(params.omit_cl_nucleotide == true && params.omit_peak_calling == true)
			"""
			wig2-to-wig.py --input ${query} --output ${query.baseName}
			sequence-extraction.py --input ${query.baseName}*.wig --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${params.percentile} --omit_cl
			"""
		else if(params.omit_cl_nucleotide == true && params.omit_peak_calling == false)
			"""
			sequence-extraction.py --input ${query} --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${params.percentile} --omit_cl
			"""
		else if(params.omit_cl_nucleotide == false && params.omit_peak_calling == true)
			"""
			wig2-to-wig.py --input ${query} --output ${query.baseName}
			sequence-extraction.py --input ${query.baseName}*.wig --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${params.percentile}
			"""
		else if(params.omit_cl_nucleotide == false && params.omit_peak_calling == false)
			"""
			sequence-extraction.py --input ${query} --reference ${reference} --output ${query.baseName}.extracted-sequences.fasta --length ${params.seq_len} --percentile ${params.percentile}
			"""
	}
}

// if sequence_extraction is performed in FASTA format, we can compute motifs
if (params.omit_sequence_extraction == false && (2*params.seq_len)+1 >= params.min_motif_width) {
    process motif_search {
    	tag {fasta.simpleName}
	    publishDir "${params.output}/motif_search/", mode: 'copy'

	    input:
	    file fasta from extracted_sequences

		output:
	    file "${fasta.baseName}_motif" optional true

		script:
		
	    """
		if [[ \$(wc -l ${fasta} | cut -f1 -d' ') -ge 4 ]]; then
			streme --oc ${fasta.baseName}_motif --p ${fasta} --dna --seed 0 --nmotifs ${params.max_motif_num} --minw ${params.min_motif_width} --maxw ${params.max_motif_width}
		fi
	    """
    }
}


if (params.omit_peak_distance == false) {
	process calculate_peak_distance {
		tag {query.simpleName}

		publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.tsv"

		input:
		file query from collected_wig_2_to_peak_distance

		output:
		file "${query.baseName}.peak-distance.tsv" into tsv_peak_distance_calculation_to_plot,tsv_peak_distance_into_output

		"""
		wig2-to-wig.py --input ${query} --output ${query.baseName}
		peak-distance.py --input ${query.baseName}_*.wig --output ${query.baseName}.peak-distance.tsv --percentile ${params.percentile} --distance ${params.distance}
		"""
	}

	process plot_peak_distance {

		publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.simpleName}.peak-distance{_full,}.png"

		input:
		file query from tsv_peak_distance_calculation_to_plot

		output:
		file "${query.simpleName}.peak-distance{_full,}.png" into png_peak_distance_into_output

		"""
		plot-distances.R --input ${query} --output ${query.simpleName}.peak-distance --type png
		"""
	}
}

process generate_barcode_barplot {
	publishDir "${params.output}/statistics", mode: 'copy'

	input:
	file(query) from log_split_experimental_barcode.first()

	output:
	file("${query.baseName}.png") into output_experimental_barcode_distribution_png

	"""
	plot_experimental_barcode_distribution.R --logs ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
	"""
}

process collect_experiments_without_alignments {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'experiments-without-alignments.txt'

	input:
	file(query) from log_experiments_without_alignments.flatten().toList()

	output:
	file("experiments-without-alignments.txt") optional true into output_log_exp_without_alignment

	"""
	if [[ ! \$(cat ${query} | wc -l) == 0 ]]; then
		cat ${query} > experiments-without-alignments.txt
	fi
	"""
}

if (params.map_to_transcripts == true){
	top_transcripts_to_collect.set{ reference_to_output }
} else {
	reference_to_collect.set{ reference_to_output }
}

process output_reference {
	publishDir "${params.output}", mode: 'copy', pattern: "${query}"

	input:
	file(query) from reference_to_output

	output:
	file(query)

	"""
	"""
}

if (params.merge_replicates == true){
	track_path_dir = Channel.value('cross-link-sites-merged/bigWig')
	big_wig_to_igv_session.filter{it[0] == 'cross-link-sites-merged'}.map{it[1]}.set{filtered_big_wig}
} else {
	track_path_dir = Channel.value('cross-link-sites/bigWig')
	big_wig_to_igv_session.map{it[1]}.set{filtered_big_wig}
}

//big_wig_to_igv_session.flatten().toList().combine(track_path_dir).set{collect_peaks}
if (params.annotation != 'NO_FILE'){
	process prepare_annotation_for_igv {
		publishDir "${params.output}", mode: 'copy', pattern: "${annotation.baseName}.sorted.gff.gz*"

		input:
		file(annotation) from annotation_to_prepare_for_igv

		output:
		val("${annotation.baseName}.sorted.gff.gz") into annotation_name_to_igv_session
		file("${annotation.baseName}.sorted.gff.gz*")

		"""
		~/software/gff3sort/gff3sort.pl ${annotation} > ${annotation.baseName}.sorted.gff
		bgzip ${annotation.baseName}.sorted.gff
		tabix ${annotation.baseName}.sorted.gff.gz
		"""
	}
} else {
	annotation_name_to_igv_session = Channel.from('NO_FILE')
}

process generate_igv_session {
	publishDir "${params.output}", mode: 'copy', pattern: 'igv-session.xml'

	input:
	file(tracks) from filtered_big_wig.flatten().toList()
	file(bam) from bam_sorted_to_igv_session.flatten().toList()
	val(track_path) from track_path_dir
	file(ref) from reference_to_igv
	val(annotation) from annotation_name_to_igv_session

	output:
	file('igv-session.xml')

	script:
	if(params.annotation == 'NO_FILE')
		"""
		generate-igv-session.py --reference ${ref} --input_path ${track_path} --tracks ${tracks} ${bam} --output igv-session.xml
		"""
	else
		"""
		generate-igv-session.py --reference ${ref} --annotation ${annotation} --input_path ${track_path} --tracks ${tracks} ${bam} --output igv-session.xml
		"""
}

process multiqc{
	publishDir "${params.output}/statistics", mode: 'move'

	input:
	file adapter from collect_statistics_adapter.first().flatten().toList()
	file qual from collect_statistics_quality_filter.first().flatten().toList()
	file qc1 from collect_statistics_qc1.first().flatten().toList()
	file qc2 from collect_statistics_qc2.first().flatten().toList()
	file split from collect_statistics_split.first().flatten().toList()
	file mapping from collect_statistics_mapping.first().flatten().toList()
	file deduplication from log_deduplicate_to_collect_statistics.first().flatten().toList()
	//file(params) from params_peak_calling_to_collect_statistics.flatten().toList().first

	output:
	file "multiqc_*" into multiqc_to_output


	"""
	multiqc .
	"""
}

workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}
