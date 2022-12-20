#!/bin/bash

nextflow run PARANOiD.nf    \
--reads    ./test_data/iCLIP_RVFV_R1.fastq  `# --reads reads_file.fastq or "reads_{1,2}.fastq" or "*.fastq"`\
--reference ./test_data/reference_RVFV.fasta    `# --reference reference_file.fasta`\
--barcodes ./test_data/barcodes-RVFV.tsv \
-profile singularity \
--config ./nextflow.config \
--resume \ 

# -with-dag flowchart.png 