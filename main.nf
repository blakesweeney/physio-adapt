#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 * Temperature Adaptation Pipeline
 *
 * This pipeline processes BacDive CSV data to:
 * 1. Convert CSV to JSONL format
 * 2. Extract GCA accession numbers
 * 3. Download genomes using NCBI datasets
 * 4. Run cmscan with Rfam covariance models
 */


/*
 * Process 1: Convert BacDive CSV to JSONL
 */
process CSV_TO_JSONL {
    tag "bacdive_conversion"
    cpus 1
    memory '2 GB'

    input:
    path csv_file

    output:
    path "bacdrive.jsonl"

    script:
    """
    bacdrive2jsonl.py ${csv_file} bacdrive.jsonl
    """
}

/*
 * Process 2: Extract GCA accession numbers from JSONL
 */
process EXTRACT_GCA_ACCESSIONS {
    tag "extract_accessions"
    cpus 1
    memory '1 GB'

    input:
    path jsonl_file

    output:
    path "gca_accessions.txt"

    script:
    """
    extract_gca_accessions.py ${jsonl_file} gca_accessions.txt --limit ${params.max_genomes}
    """
}

/*
 * Process 3: Download Rfam covariance models
 */
process DOWNLOAD_RFAM {
    tag "download_rfam"
    cpus 1
    memory '2 GB'
    time '1 h'

    output:
    tuple path("Rfam.cm*"), path("Rfam.clanin")

    script:
    """
    # Download and decompress Rfam.cm
    curl -L ${params.rfam_cm_url} -o Rfam.cm.gz
    gunzip Rfam.cm.gz

    # Download Rfam.clanin
    curl -L ${params.rfam_clanin_url} -o Rfam.clanin

    # Press the Rfam.cm file for faster searching (creates index files)
    cmpress Rfam.cm
    """
}

/*
 * Process 4: Download genome sequences using NCBI datasets
 */
process DOWNLOAD_GENOME {
    tag "${accession}"
    cpus 1
    memory '2 GB'
    time '1 h'
    errorStrategy 'ignore'

    input:
    val accession

    output:
    tuple val(accession), path("${accession}.fa"), optional: true

    script:
    """
    # Download genome using datasets command
    datasets download genome accession ${accession} \
        --include genome \
        --filename ${accession}.zip

    # Extract the genome sequence
    unzip -q ${accession}.zip

    # Find and concatenate all genomic fasta files
    find ncbi_dataset/data -name "${accession}*_genomic.fna" | xargs -I {} cp {} ${accession}.fa

    # Check if we got a valid genome file
    if [ ! -s ${accession}.fa ]; then
        echo "Warning: No genome sequence found for ${accession}" >&2
        exit 1
    fi
    """
}

/*
 * Process 5: Run cmscan with Rfam models on each genome
 */
process CMSCAN {
    tag "${accession}"
    cpus 4
    memory '8 GB'
    time '4 h'
    publishDir "${params.outdir}/cmscan", mode: 'copy'

    input:
    tuple val(accession), path(genome_fna), path(rfam_cm), path(clanin)

    output:
    path "${accession}.tblout"

    script:
    """
    # Run cmscan with Rfam cutoffs and clan information
    cmscan --cut_ga \
        --rfam \
        --nohmmonly \
        --tblout ${accession}.tblout \
        --fmt 2 \
        --clanin ${clanin} \
        --cpu ${task.cpus} \
        Rfam.cm \
        ${genome_fna}
    """
}

/*
 * Workflow
 */
workflow {
    // Convert CSV to JSONL
    Channel.fromPath(params.input_csv, checkIfExists: true) | set { bacdrive }

    DOWNLOAD_RFAM | set { rfam }

    bacdrive \
    | CSV_TO_JSONL \
    | EXTRACT_GCA_ACCESSIONS \
    | splitText \
    | map { it.trim() } \
    | filter { it.length() > 0 } \
    | DOWNLOAD_GENOME \
    | set { genomes }

    genomes \
    | combine(rfam) \
    | CMSCAN
}
