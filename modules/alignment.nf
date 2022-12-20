// processes for alignment 

process build_index_bowtie {

	input:
	path ref // from reference_to_mapping

	output:
	tuple file("${ref}"), file("${ref}.*") // into bowtie_index_build_to_mapping

	"""
	bowtie2-build ${ref} ${ref}
	"""
}

process mapping_bowtie{
	tag {query.simpleName}

	input:
	tuple file(ref), file(index) // from bowtie_index_build_to_mapping.first()
	path query // from fastq_merge_preprocessed_to_alignment

	output:
	path "${query.baseName}.bam" ,emit: bam_mapping // into (bam_mapping_to_filter_empty, bam_mapping_to_sort)
	path "${query.simpleName}.statistics.txt", emit: collect_statistics_mapping   // into collect_statistics_mapping

	"""
	bowtie2 --no-unal -q -p ${task.cpus} -U ${query} -x ${ref} 2> ${query.simpleName}.statistics.txt | samtools view -bS - > ${query.baseName}.bam
	"""
}

	//Version 2.7.3a
	process build_index_STAR {

	input:
	path referenceGenome // from reference_to_mapping
	path gtf // from sjdbGTFfile

	output:
	path index // into star_index_build_to_mapping

	script: 
	if(params.annotation == 'NO_FILE')
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} 
		"""
	else
		"""
		mkdir index
		STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ./index --genomeFastaFiles ${referenceGenome} --sjdbGTFpath ${gtf}
		"""
}


process mapping_STAR{
	tag {query.simpleName}

	input:
	tuple file(query), file(indexDir) // from fastq_merge_preprocessed_to_alignment.combine(star_index_build_to_mapping)

	output:
	file("${query.baseName}.Aligned.sortedByCoord.out.bam") , emit bam_mapping // into bam_mapping_to_filter_empty
	file("${query.baseName}.Log.*")  ,emit: collect_statistics_mapping //  into collect_statistics_mapping

	"""
	STAR --runThreadN ${task.cpus} --genomeDir ${indexDir} --readFilesIn ${query} --outFileNamePrefix ${query.baseName}. --alignEndsType Extend5pOfRead1 --outSAMtype BAM SortedByCoordinate
	"""
}

