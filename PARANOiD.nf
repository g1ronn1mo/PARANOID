#!/usr/bin/env nextflow


nextflow.enable.dsl=2		// Enable nextflow DSL version 2

import java.nio.file.*
import groovy.io.FileType

// import processes for preprocessing modules 
include {
	quality_control
	quality_control_2
	adapter_removal
	quality_filter
} from './modules/preprocessing.nf'


// import processes for barcode wrangling
include {
	extract_rnd_barcode
	check_barcode_file
	split_exp_barcode
	get_length_exp_barcode
	remove_exp_barcode
	merge_preprocessed_reads
	generate_barcode_barplot
} from './modules/barcode.nf'


// import processes for alignment 
include {
	build_index_bowtie
	mapping_bowtie
	build_index_STAR
	mapping_STAR
} from './modules/alignment.nf'


process filter_empty_bams{
	tag {query.simpleName}

	input:
	path query // from bam_mapping_to_filter_empty

	output:
	path "${query.simpleName}.no_alignments.txt", emit: log_experiments, optional:  true // into log_experiments_without_alignments
	path "${query.baseName}.filtered.bam", emit: bam_filter_empty, optional:  true // into bam_filter_empty_to_split

	"""
	if [[ \$(samtools view ${query} | wc -l) == 0 ]]; then
		echo "${query.simpleName} contains no alignable reads" > ${query.simpleName}.no_alignments.txt
	else
		mv ${query} ${query.baseName}.filtered.bam
	fi
	"""
}


process bam_filter_empty{
	tag {query.simpleName}

	input:
	file(query) // from bam_filter_empty_to_split

	output:
	file("*.bam") // into bam_split_to_sort

	"""
	bamtools split -in ${query} -reference
	"""
}

process sort_bam{
	tag {query.simpleName}

	input:
	path query // from bam_split_to_sort.flatten()

	output:
	tuple file("${query.baseName}.sorted.bam"), file("${query.baseName}.sorted.bam.bai")  // into bam_sort_to_deduplicate

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
	tuple file(query), file(index) // from bam_sort_to_deduplicate

	output:
	path "${query.baseName}.deduplicated.bam" ,emit: bam_deduplicated // into (bam_deduplicate_to_sort ,bam_deduplicate_to_index)
	path "${query.baseName}.deduplicated.log*" ,emit: log_deduplicate // into (log_deduplicate_to_collect_statistics,log_deduplicate_to_output)
	
	"""
	umi_tools dedup --random-seed=42 -I ${query} --output-stats ${query.baseName}.deduplicated.log -S ${query.baseName}.deduplicated.bam
	"""
}

process merge_deduplicated_bam {
	tag {name}
	input:
	tuple val(name), file(query) // from bam_dedup_sort_to_merge

	output:
	path ("${name}.bam"), emit: bam_merge // into (bam_merge_to_calculate_crosslinks, bam_merge_to_extract_transcripts, bam_merge_to_pureCLIP) 
	// path ("${name}.bam"), emit: bam_alignment  // into bam_alignment_to_sort

	"""
	samtools merge -c -p ${name}.bam ${query}
	"""
}

process sort_and_index_alignment{
	tag {query.simpleName}
	publishDir "${params.output}/alignments", mode: 'copy', pattern: "${query.simpleName}.sorted.bam*"

	input:
	path (query) // from bam_alignment_to_sort

	output:
	path ("${query.simpleName}.sorted.bam") , emit: bam_sorted_to_igv_session   // into bam_sorted_to_igv_session
	path ("${query.simpleName}.sorted.bam*") 

	"""
	samtools sort ${query} > ${query.simpleName}.sorted.bam
	samtools index ${query.simpleName}.sorted.bam
	"""
}


process count_hits {
	tag {bam.simpleName}

	publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

	input:
	path bam // from bam_merge_to_extract_transcripts

	output:
	path ("${bam.baseName}.hits.tsv") // into tsv_count_to_get_hits

	"""
	samtools view ${bam} | cut -f3 | sort | uniq -c | sort -nr > ${bam.baseName}.hits.tsv
	"""
}

	process get_top_hits {

		publishDir "${params.output}/transcripts/overview-hits", mode: 'copy'

		input:
		path tsv // from tsv_count_to_get_hits.flatten().toList()

		output:
		path ("transcript-targets-top${params.number_top_transcripts}.txt") // into (txt_get_hits_to_extract_alignments,txt_get_hits_to_extract_sequences)

		"""
		head -${params.number_top_transcripts} -q *.tsv  | rev | cut -f1 -d' ' | rev | sort | uniq > transcript-targets-top${params.number_top_transcripts}.txt
		"""
	}

	process bam_deduplicate {
		tag {bam.simpleName}

		input:
		path bam // from bam_deduplicate_to_index

		output:
		tuple file("${bam.baseName}.sorted.bam"), file("${bam.baseName}.sorted.bam.bai") // into bam_index_to_extract_alignments

		"""
		samtools sort ${bam} -o ${bam.baseName}.sorted.bam
		samtools index ${bam.baseName}.sorted.bam
		"""
	}

	process extract_top_alignments {
		tag {bam.simpleName}

		input:
		tuple file(bam), file(bai), file(txt_alignments) // from bam_index_to_extract_alignments.combine(txt_get_hits_to_extract_alignments)

		output:
		path ("${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam") // into remove_newlines_to_calc_crosslink

		"""
		samtools view -hb ${bam} `cat ${txt_alignments}` > ${bam.simpleName}_filtered_top${params.number_top_transcripts}.bam
		"""
	}

	process remove_newlines {

		input:
		path ref // from reference_to_extract_transcripts

		output:
		path ("${ref.baseName}.removed_newlines.fna") // into fasta_rm_newline_to_extract_sequences

		"""
		awk '!/^>/ {printf "%s", \$0; n = "\\n" } /^>/ { print n \$0; n = "" } END { printf "%s", n }' ${ref} > ${ref.baseName}.removed_newlines.fna
		"""
	}

	process extract_top_transcript_sequences {
		tag {txt_sequences.simpleName}

		publishDir "${params.output}/transcripts", mode: 'copy'

		input:
		tuple path (txt_sequences),path (ref) // from txt_get_hits_to_extract_sequences.combine(fasta_rm_newline_to_extract_sequences)

		output:
		path ("${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna") // into (top_transcripts_to_collect,top_transcripts_to_strand_preference)


		"""
		egrep -A1 --no-group-separator -f ${txt_sequences} ${ref} > ${ref.simpleName}.top${params.number_top_transcripts}_transcripts.fna
		"""
	}


process sort_bam_before_strand_pref {
	tag {query.baseName}

	input:
	path (query) // from collected_bam_files_to_sort

	output:
	path ("${query.baseName}.sorted.bam") // into bam_sort_to_group

	"""
	samtools sort ${query} > ${query.baseName}.sorted.bam
	"""
}


process determine_strand_preference {
	tag {name}
	publishDir "${params.output}/strand-distribution", mode: 'copy', pattern: "${name}.strand_proportion.txt"

	input:
	tuple val(name),file(query),file(reference) // from grouped_bam_to_strand_preference.combine(choose_correct_reference_file)

	output:
	path ("${name}.strand_proportion.txt") // into txt_determine_strand_to_visualize

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
	path (strand) // from txt_determine_strand_to_visualize

	output:
	path ("${strand.simpleName}.png")

	"""
	visualize_strand_distribution.R --input ${strand} --output ${strand.simpleName} --type png
	"""
}

process get_chromosome_sizes{
	input:
	path (ref) // from reference_to_chrom_sizes

	output:
	path ("${ref.simpleName}.chromosome_sizes.txt") // into (chrom_sizes_to_cross_link_calculation,chrom_sizes_to_bigWig)

	"""
	samtools faidx ${ref}
	cut -f1,2 ${ref}.fai > ${ref.simpleName}.chromosome_sizes.txt
	"""
}

process calculate_crosslink_sites{
	tag {query.simpleName}
	publishDir "${params.output}/cross-link-sites/wig", mode: 'copy', pattern: "${query.simpleName}_{forward,reverse}.wig"

	input:
	tuple file(query), file(chrom_sizes) // from collected_bam_files.combine(chrom_sizes_to_cross_link_calculation)

	output:
	path "${query.simpleName}.wig2", emit: wig_calculate_crosslink,   optional: true // into (wig_calculate_crosslink_to_group_samples,wig2_calc_cl_sites_to_split)
	tuple val("cross-link-sites"), file("${query.simpleName}_forward.wig"), file("${query.simpleName}_reverse.wig"), emit: wig_cross_link_site, optional: true // into wig_cross_link_sites_to_transform
	path "${query.simpleName}_{forward,reverse}.wig", emit:  wig_calculate_cl_sites , optional: true // into wig_calculate_cl_sites_to_split_strand

	"""
	create-wig-from-bam.py --input ${query} --mapq ${params.mapq} --chrom_sizes ${chrom_sizes} --output ${query.simpleName}.wig2
	if [[ -f "${query.simpleName}.wig2" ]]; then
		wig2-to-wig.py --input ${query.simpleName}.wig2 --output ${query.simpleName}
	fi
	"""
}

//// from here on only further analyses
process split_wig2_for_correlation{
	tag {query.simpleName}

	input:
	file(query) // from wig2_calc_cl_sites_to_split

	output:
	file("${query.simpleName}_forward.wig"), emit: split_wig_forward optional true // into split_wig_forward_to_correlation
	file("${query.simpleName}_reverse.wig"), emit:  split_wig_reverse optional true // into split_wig_reverse_to_correlation

	"""
	wig2-to-wig.py --input ${query} --output ${query.simpleName}
	"""

}

/*process calc_wig_correlation{
	tag {name}
	echo true
	cache false

	input:
	tuple val(name),file(query) // from group_forward

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
	tuple name, file(query) // from wig2_grouped_samples_to_merge

	output:
	path "${name}.wig2" // into collected_wig_files
	tuple val("cross-link-sites-merged"), file("${name}_forward.wig"), file("${name}_reverse.wig"), emit: wig_merged_cross_link_sites // into wig_merged_cross_link_sites_to_transform
	path "${name}_{forward,reverse}.wig" // into output

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



process wig_to_bigWig{
	tag {forward.simpleName}
	publishDir "${params.output}/${out_dir}/bigWig", mode: 'copy', pattern: "*.bw"

	input:
	tuple val(out_dir), file(forward), file(reverse), file(chrom_sizes) // from wig_cross_link_sites_to_bigWig.combine(chrom_sizes_to_bigWig)

	output:
	file("*.bw") optional true
	tuple val(out_dir), file("*.bw"), emit: big_wig, optional: true// into big_wig_to_igv_session
	tuple val(out_dir), file("${reverse.baseName}.bw"), emit: big_wig_reverse, optional: true // into big_wig_reverse_to_convert_to_bedgraph
	tuple val(out_dir), file("${forward.baseName}.bw"), emit: big_wig_forward, optional: true // into big_wig_forward_to_convert_to_bedgraph

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
	tuple val(out_dir), file(bigWig) // from big_wig_forward_to_convert_to_bedgraph.mix(big_wig_reverse_to_convert_to_bedgraph)

	output:
	file("*.bedgraph")

	"""
	bigWigToBedGraph ${bigWig} ${bigWig.baseName}.bedgraph
	"""
}


process index_for_peak_calling {
	tag{query.simpleName}

	input:
	path (query) // from bam_merge_to_pureCLIP

	output:
	tuple path ("${query.simpleName}.sorted.bam"), path ("${query.simpleName}.sorted.bam.bai") // into bambai_index_to_peak_calling
	
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
	tuple file(bam), file(bai), file(ref) // from bambai_index_to_peak_calling.combine(reference_to_pureCLIP)

	output:
	path ("${bam.simpleName}.pureCLIP_crosslink_{sites,regions}.bed")
	path ("${bam.simpleName}.pureCLIP_crosslink_sites.bed"), emit:bed_peak_calling // into bed_peak_calling_to_sequence_extraction
	path ("${bam.simpleName}.pureCLIP_crosslink_sites.params"), emit:params_peak_calling // into params_peak_calling_to_collect_statistics

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


process split_wig_2_for_peak_height_hist {
	tag {query.simpleName}

	input:
	path (query) // from collect_wig_2_to_peak_height_histogram

	output:
	tuple val("${query.simpleName}"), file("${query.simpleName}_forward.wig"), file("${query.simpleName}_reverse.wig") // into split_wig2_to_generate_peak_height_histogram

	"""
	wig2-to-wig.py --input ${query} --output ${query.simpleName}
	"""
}

process generate_peak_height_histogram {
	tag {query}
	publishDir "${params.output}/peak_height_distribution", mode: 'copy'

	input:
	tuple val(query), file(forward), file(reverse) // from split_wig2_to_generate_peak_height_histogram

	output:
	path("${query}.png")

	"""
	generate-peak-height-histogram.R --input . --output ${query} --type png --color "${params.color_barplot}" --percentile ${params.percentile}
	"""
}


process wig_to_bam {
	tag {query.simpleName}

	input:
	path(query) // from collected_wig_2_to_RNA_subtypes_distribution

	output:
	path("${query.baseName}.bam") // into bam_convert_to_feature_counts

	"""
	wig-to-bam.py --input ${query} --output ${query.baseName}.bam
	"""
}

process feature_counts {
	tag {query.simpleName}

	input:
	tuple file(query), val(rna_subtypes), file(annotation) // from bam_convert_to_feature_counts.combine(rna_subtypes_to_feature_counts).combine(annotation_to_RNA_subtypes_distribution)

	output:
	path("${query.simpleName}.${rna_subtypes}.tsv") // into tsv_feature_counts_to_sort


	"""
	featureCounts -T ${task.cpus} -t ${rna_subtypes} -g ${params.gene_id} -a ${annotation} -R CORE -M -o ${query.simpleName} ${query}
	mv ${query}.featureCounts ${query.simpleName}.${rna_subtypes}.tsv
	"""
}


process get_RNA_subtypes_distribution {
	tag {name}

	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	val(subtypes) // from rna_subtypes_to_distribution.collect()
	tuple val(name), file(query) // from tsv_sort_to_calculate_distribution

	output:
	path("${name}.subtypes_distribution.tsv"), emit: output_rna_subtypes_tsv  // into (output_rna_subtypes_tsv,tsv_rna_subtypes_distribution_to_barplot)
	path("${name}.subtype.log"), emit:log_subtype_warnings  ,optional: true // into log_subtype_warnings

	script:
	subtypes_as_string = subtypes.join(' ')
	"""
	calc-RNA-subtypes-distribution.py --input ${query} --rna_subtypes ${subtypes_as_string} --output ${name}.subtypes_distribution.tsv > ${name}.subtype.log
	"""
}

process generate_RNA_subtypes_barplot {
	publishDir "${params.output}/RNA_subtypes", mode: 'copy'

	input:
	path(query) // from tsv_rna_subtypes_distribution_to_barplot

	output:
	path("${query.baseName}.png") // into output_rna_subtypes_png

	"""
	RNA_subtypes_barcharts.R --input ${query} --output ${query.baseName} --type png --color "${params.color_barplot}"
	"""
}

process collect_subtype_analysis_errors {
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'subtype-analysis-warnings.txt'

	input:
	path(query) // from log_subtype_warnings.flatten().toList()

	output:
	file("subtype-analysis-warnings.txt") optional true // into output_log_subtype_analysis_warnings

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


process sequence_extraction {
	tag {query.simpleName}

	publishDir "${params.output}/extracted_sequences", mode: 'copy', pattern: "*.extracted-sequences.*"

	input:
	tuple file(query),file(reference) // from peaks_for_sequence_extraction.combine(reference_to_extract_sequences)

	output:
	path "*.extracted-sequences.fasta", emit:extracted_sequences, optional: true // into extracted_sequences
	path "*.extracted-sequences.*",emit:extracted_sequences_output , optional: true // into extracted_sequences_output
	
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


// if sequence_extraction is performed in FASTA format, we can compute motifs
process motif_search {
    	tag {fasta.simpleName}
	    publishDir "${params.output}/motif_search/", mode: 'copy'

	    input:
	    path fasta // from extracted_sequences

		output:
	    path "${fasta.baseName}_motif", optional: true

		script:
		
	    """
		if [[ \$(wc -l ${fasta} | cut -f1 -d' ') -ge 4 ]]; then
			streme --oc ${fasta.baseName}_motif --p ${fasta} --dna --seed 0 --nmotifs ${params.max_motif_num} --minw ${params.min_motif_width} --maxw ${params.max_motif_width}
		fi
	    """
}




process calculate_peak_distance {
	tag {query.simpleName}

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.baseName}.peak-distance.tsv"

	input:
	path query // from collected_wig_2_to_peak_distance

	output:
	path "${query.baseName}.peak-distance.tsv" // into tsv_peak_distance_calculation_to_plot,tsv_peak_distance_into_output

	"""
	wig2-to-wig.py --input ${query} --output ${query.baseName}
	peak-distance.py --input ${query.baseName}_*.wig --output ${query.baseName}.peak-distance.tsv --percentile ${params.percentile} --distance ${params.distance}
	"""
}

process plot_peak_distance {

	publishDir "${params.output}/peak_distance", mode: 'copy', pattern: "${query.simpleName}.peak-distance{_full,}.png"

	input:
	path query // from tsv_peak_distance_calculation_to_plot

	output:
	path "${query.simpleName}.peak-distance{_full,}.png" // into png_peak_distance_into_output

	"""
	plot-distances.R --input ${query} --output ${query.simpleName}.peak-distance --type png
	"""
}




process collect_experiments_without_alignments {
	tag {query.simpleName}
	publishDir "${params.output}/statistics", mode: 'copy', pattern: 'experiments-without-alignments.txt'

	input:
	path(query) // from log_experiments_without_alignments.flatten().toList()

	output:
	path("experiments-without-alignments.txt"), optional :true // into output_log_exp_without_alignment

	"""
	if [[ ! \$(cat ${query} | wc -l) == 0 ]]; then
		cat ${query} > experiments-without-alignments.txt
	fi
	"""
}

process output_reference {
	publishDir "${params.output}", mode: 'copy', pattern: "${query}"

	input:
	path(query) // from reference_to_output

	output:
	path(query)

	"""
	"""
}

//big_wig_to_igv_session.flatten().toList().combine(track_path_dir).set{collect_peaks}

process prepare_annotation_for_igv {
	publishDir "${params.output}", mode: 'copy', pattern: "${annotation.baseName}.sorted.gff.gz*"

	input:
	path(annotation) // from annotation_to_prepare_for_igv

	output:
	val("${annotation.baseName}.sorted.gff.gz") , emit: annotation_name // into annotation_name_to_igv_session
	path("${annotation.baseName}.sorted.gff.gz*")

	"""
	~/software/gff3sort/gff3sort.pl ${annotation} > ${annotation.baseName}.sorted.gff
	bgzip ${annotation.baseName}.sorted.gff
	tabix ${annotation.baseName}.sorted.gff.gz
	"""
}



process generate_igv_session {
	publishDir "${params.output}", mode: 'copy', pattern: 'igv-session.xml'

	input:
	path(tracks) // from filtered_big_wig.flatten().toList()
	path(bam) // from bam_sorted_to_igv_session.flatten().toList()
	val(track_path) // from track_path_dir
	file(ref) // from reference_to_igv
	val(annotation) // from annotation_name_to_igv_session

	output:
	path('igv-session.xml')

	script:
	if(params.annotation == 'NO_FILE')
		"""
		generate-igv-session.py --reference ${reference} --input_path ${track_path} --tracks ${tracks} ${bam} --output igv-session.xml
		"""
	else
		"""
		generate-igv-session.py --reference ${reference} --annotation ${annotation} --input_path ${track_path} --tracks ${tracks} ${bam} --output igv-session.xml
		"""
}

process multiqc{
	publishDir "${params.output}/statistics", mode: 'move'

	input:
	path adapter // from collect_statistics_adapter.first().flatten().toList()
	path qual // from collect_statistics_quality_filter.first().flatten().toList()
	path qc1 // from collect_statistics_qc1.first().flatten().toList()
	path qc2 // from collect_statistics_qc2.first().flatten().toList()
	path split // from collect_statistics_split.first().flatten().toList()
	path mapping // from collect_statistics_mapping.first().flatten().toList()
	path deduplication // from log_deduplicate_to_collect_statistics.first().flatten().toList()
	// file(params) // from params_peak_calling_to_collect_statistics.flatten().toList().first

	output:
	path "multiqc_*" // into multiqc_to_output

	"""
	multiqc .
	"""
}

process split_bam_by_chromosome{
		tag {query.simpleName}

		input:
		file(query) // from bam_filter_empty_to_split

		output:
		file("*.bam") // into bam_split_to_sort

		"""
		bamtools split -in ${query} -reference
		"""
}


sjdbGTFpath = file(params.annotation) // GTF path containing annotation information
reference = Channel.fromPath( params.reference )		//FASTA path containing reference sequence(s)
input_reads = Channel.fromPath( params.reads )			//FASTQ file(s) containing reads
rna_subtypes = Channel.of(params.rna_subtypes.split(',')) 


// ##################################################################################################
// ####################################            log             ##################################
// ##################################################################################################

log.info """
BIOCORE@CRG - N F TESTPIPE  ~  version ${params.version}
=============================================
reads                           : ${params.reads}
reference                       : ${params.reference}
"""

// ##################################################################################################
// ####################################        Workflows           ##################################
// ##################################################################################################

workflow preprocessing {
	take: input_reads

	main:
		if (params.speed) {
			input_reads.splitFastq( by: params.split_fastq_by, file:true ).set{input_reads}
		} 

		quality_control(input_reads)  // into out1, collect_statistics_qc1 
		adapter_removal(input_reads)  // emit: trimmed_fq     emit: trimming_report
		quality_filter(adapter_removal.out.trimmed_fq) // into reads_qualityFilter // into collect_statistics_quality_filter
		quality_control_2(quality_filter.out.filter_fq)  //into fastq_quality_filter_to_quality_control_2

	emit:
		// data for qc
		mqc_quality_control  = quality_control.out
		mqc_adapter_removal  = adapter_removal.out.trimming_report
		mqc_quality_control_2 = quality_control_2.out
		mqc_quality_filter   = quality_filter.out.filter_report

		// data for downstream analysis
		reads_qualityFilter = quality_filter.out.filter_fq
	}

workflow barcode {
	take: reads_qualityFilter
	main:
		//essential inputs
		val_barcode_pattern = Channel.of( params.barcode_pattern )
		barcode_file = Channel.fromPath( params.barcodes ) 	//TSV path containing experiment names and the corresponding experiemental barcode sequence
		
			extract_rnd_barcode(reads_qualityFilter)
		check_barcode_file(barcode_file)
		// extract_rnd_barcode.out.rnd_barcode.view { "extract_rnd_barcode.out.rnd_barcode.view : $it" } // rndBarcode.fastq
		// check_barcode_file.out.view { "check_barcode_file.out.view: $it" }
		split_exp_barcode( extract_rnd_barcode.out.rnd_barcode.combine(check_barcode_file.out) )
		get_length_exp_barcode(val_barcode_pattern)
		fastq_remove_barcode_to_collect = remove_exp_barcode(split_exp_barcode.out.fastq_split_barcode_to_remove_exp_barcode.flatten().filter{ it.size() > 0 }.combine(get_length_exp_barcode.out))

		remove_exp_barcode.out
		.map{file -> tuple(file.name - ~/\.[\w.]+.fastq$/,file)}
		.groupTuple()
		.set{ fastq_collect_preprocessed_to_merge }
		merge_preprocessed_reads(fastq_collect_preprocessed_to_merge)

		output_experimental_barcode_distribution_png = generate_barcode_barplot(split_exp_barcode.out.log.first())

	emit:
		// data for qc report
		split_exp_barcode_qc = split_exp_barcode.out.log

		// data for downstream analysis
		fastq_merge_preprocessed = merge_preprocessed_reads.out
	}

workflow alignment {
	//essential inputs    
	take: fastq_merge_preprocessed

	main:
		if ( params.domain == 'pro' || params.map_to_transcripts == true)
		{
		//bowtie2 version 2.3.5.1
		bowtie_index_build = build_index_bowtie(reference)
		mapping_bowtie(bowtie_index_build.first(), fastq_merge_preprocessed)
			mapping = mapping_bowtie.out.bam_mapping
			stat = mapping_bowtie.out.collect_statistics_mapping
		} 
		else if ( params.domain == 'eu' ) 
		{
		// 	STAR
		star_index_build_to_mapping = build_index_STAR( reference ,  sjdbGTFfile)
		mapping_STAR (fastq_merge_preprocessed.combine(star_index_build))
			mapping = mapping_STAR.out.bam_mapping
			stat = mapping_STAR.out.collect_statistics_mapping
		}
		
		emit:
			mapping
			stat
}

	
// Todo: Split the process_bam Workflow into two workflows: one for the processing of the bam files and one for the processing of the bigwig files
workflow process_bam {
	take: mapping_in
	main:
		filter_empty_bams(mapping_in)
		collect_experiments_without_alignments(filter_empty_bams.out.log_experiments.flatten().toList())
		if( params.map_to_transcripts == false && params.speed ) {
			bam_split = split_bam_by_chromosome(mapping_in)
		} else   {
		bam_split = mapping_in
		}

		bam_sort = sort_bam (bam_split.flatten()) 

		bam_deduplicated = 		deduplicate(bam_sort).bam_deduplicated 
		log_deduplicate = deduplicate.out.log_deduplicate

		bam_deduplicated
		.map{file -> tuple(file.name - ~/\.[\w.]+.bam$/,file)}
		.groupTuple()
		.set{ bam_dedup_sort }

		bam_merge_deduplicated = merge_deduplicated_bam(bam_dedup_sort) 
				
		if (params.map_to_transcripts == true){ 
			tsv_count = count_hits(bam_merge_deduplicated.bam_merge)
			txt_get_hits = get_top_hits( tsv_count.flatten().toList())
			index_alignments
			bam_extract_alignments = extract_top_alignments(bam_deduplicated.combine(txt_get_hits))
			fasta_rm_newline = remove_newlines(reference)
			top_transcripts = extract_top_transcript_sequences(txt_get_hits.combine(fasta_rm_newline))
			collected_bam_files = bam_extract_alignments

			top_transcripts.set{choose_correct_reference_file}

		} else {
			collected_bam_files  = bam_merge_deduplicated.bam_merge
			reference.set{choose_correct_reference_file}
		}

		if(params.merge_replicates == true){
			bam_sort
			.map{file -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*(.sorted)?.bam$/,file)} 
			.groupTuple()
			.set{grouped_bam}
		} else {
			bam_sort
			.map{file -> tuple(file.name - ~/(_filtered_top\d*)?(.sorted)?.bam$/,file)}
			.set{grouped_bam}
		} 

		// Todo: generate_igv_session does not work!
		//  determine_strand_preference(grouped_bam.combine(choose_correct_reference_file)) | visualize_strand_preference

		chrom_sizes = get_chromosome_sizes(reference)
		calculate_crosslink_sites(collected_bam_files.combine(chrom_sizes))

		if( params.merge_replicates == true ){
			// according to their experiment
			calculate_crosslink_sites.wig_calculate_crosslink
			.map{path -> tuple(file.name - ~/_rep_\d*(_filtered_top)?\d*.wig2$/,file)} 
			.groupTuple()
			.set{wig2_grouped_samples}

			merge_wigs(wig2_grouped_samples)

			split_wig2_for_correlation ( calculate_crosslink_sites.wig_calculate_crosslink)

			split_wig2_for_correlation.split_wig_forward
			.map{path -> tuple(file.name - ~/_rep_\d*_forward.wig$/,file)} 
			.groupTuple()
			.set{group_forward}

			// Transformation of cross-link sites to different formats
			calculate_crosslink_sites.wig_cross_link_site.mix( merge_wigs.wig_merged_cross_link_sites)
			.set{wig_cross_link_sites}

		} else {

			calculate_crosslink_sites.out.wig_calculate_crosslink
			.set{collected_wig_files}
			
			calculate_crosslink_sites.out.wig_cross_link_site
			.set{wig_cross_link_sites}
		}

		wig_to_bigWig( wig_cross_link_sites.combine(chrom_sizes))
		bigWig_to_bedgraph( wig_to_bigWig.out.big_wig_forward.mix(wig_to_bigWig.out.big_wig_reverse)) 
		bam_sort_to_group = sort_bam_before_strand_pref(collected_bam_files)



		if (params.omit_peak_calling == false){
			bambai_index = index_for_peak_calling(bam_merge_deduplicated.bam_merge)
			pureCLIP(bambai_index.combine(reference))
		}

		split_wig_2_for_peak_height_hist(collected_wig_files) | generate_peak_height_histogram 


		if ( params.annotation != 'NO_FILE'){
			wig_to_bam(collected_wig_2) | feature_counts.set{tsv_feature_counts}
			
			tsv_feature_counts
			.map{file-> tuple(file.name - ~/\.[\w.]+.tsv$/,file)}
			.groupTuple()
			.set{tsv_sort}

			get_RNA_subtypes_distribution(rna_subtypes, tsv_sort )
			generate_RNA_subtypes_barplot (get_RNA_subtypes_distribution.out.output_rna_subtypes_tsv )
			collect_subtype_analysis_errors(get_RNA_subtypes_distribution.log_subtype_warnings.flatten().toList())
		}


		if (params.omit_sequence_extraction == false){
			if(params.omit_peak_calling == false){
				pureCLIP.out.bed_peak_calling.set{peaks_for_sequence_extraction}
			} else {
				collected_wig_2.set{peaks_for_sequence_extraction}
			}
			sequence_extraction(peaks_for_sequence_extraction.combine(reference)) 
		}  
		
		if (params.omit_sequence_extraction == false && (2*params.seq_len)+1 >= params.min_motif_width) {
			motif_search(sequence_extraction.out.extracted_sequences)
		}

		if (params.omit_peak_distance == false) {
			calculate_peak_distance(collected_wig_files) | plot_peak_distance
		}

		if (params.merge_replicates == true){
			track_path_dir = Channel.value('cross-link-sites-merged/bigWig')
			wig_to_bigWig.out.big_wig.filter{it[0] == 'cross-link-sites-merged'}.map{it[1]}.set{filtered_big_wig}
		} else {
			track_path_dir = Channel.value('cross-link-sites/bigWig')
			wig_to_bigWig.out.big_wig.map{it[1]}.set{filtered_big_wig}
		} 

		if (params.annotation != 'NO_FILE'){
			annotation_file = Channel.fromPath( params.annotation )
			prepare_annotation_for_igv(annotation_file).out.annotation_name
		} else {
		annotation_name= Channel.of('NO_FILE')
		}
		sort_and_index_alignment(bam_merge_deduplicated.bam_merge)
		
		// Todo: generate_igv_session does not work!
		// generate_igv_session(
		// 	filtered_big_wig.flatten().toList() , 
		// 	sort_and_index_alignment.out.bam_sorted_to_igv_session.flatten().toList(), 
		// 	track_path_dir, 
		// 	reference, 
		// 	annotation_name
		// 	) 

	emit:
		// data for qc report
		log_deduplicate 
}

// ####################################     Workflow starts here   ##################################

workflow {
	
	preprocessing(input_reads)
	barcode(preprocessing.out.reads_qualityFilter)
	alignment(barcode.out.fastq_merge_preprocessed)	
	process_bam(alignment.out.mapping)
	output_reference(reference)


	multiqc( 
		preprocessing.out.mqc_adapter_removal.first().flatten().toList(),  // from collect_statistics_adapter.first().flatten().toList() x
		preprocessing.out.mqc_quality_filter.first().flatten().toList(),  // from collect_statistics_quality_filter.first().flatten().toList() x 
		preprocessing.out.mqc_quality_control.first().flatten().toList(), // from collect_statistics_qc1.first().flatten().toList()
		preprocessing.out.mqc_quality_control_2.first().flatten().toList(), // from collect_statistics_qc2.first().flatten().toList()
		alignment.out.mapping.first().flatten().toList(),  // from collect_statistics_split.first().flatten().toList()
		barcode.out.split_exp_barcode_qc.first().flatten().toList(),  // from collect_statistics_mapping.first().flatten().toList()
		process_bam.out.log_deduplicate.first().flatten().toList() // file(params) // from params_peak_calling_to_collect_statistics.flatten().toList().first
	)
}


workflow.onComplete{
	println "Pipeline completed at: $workflow.complete"
	println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
	println "Something went wrong :("
	println "Pipeline execution stopped with following error message: ${workflow.errorMessage}"
}