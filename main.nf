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
    publishDir "${params.outdir}/bacdrive", mode: 'copy'

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
    extract_gca_accessions.py ${jsonl_file} gca_accessions.txt
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
 * Process 3B: Download Rfam metadata files for database building
 */
process DOWNLOAD_RFAM_METADATA {
    tag "download_rfam_metadata"
    cpus 1
    memory '2 GB'
    time '1 h'
    publishDir "${params.outdir}/rfam_metadata", mode: 'copy'

    output:
    tuple path("family.txt.gz"), path("clan.txt.gz"), path("clan_membership.txt.gz")

    script:
    """
    # Download Rfam database files from FTP
    curl -L https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/database_files/family.txt.gz -o family.txt.gz
    curl -L https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/database_files/clan.txt.gz -o clan.txt.gz
    curl -L https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/database_files/clan_membership.txt.gz -o clan_membership.txt.gz
    """
}

/*
 * Process 3C: Extract Rfam models by RNA type for below-threshold scanning
 */
process EXTRACT_MODELS_BY_TYPE {
    tag "extract_models"
    cpus 2
    memory '8 GB'
    time '2 h'

    input:
    tuple path(rfam_cm), path(clanin)
    tuple path(family_txt_gz), path(clan_txt_gz), path(clan_membership_txt_gz)
    path rna_types_file

    output:
    tuple path("subset.cm*"), path("selected_families.txt")

    script:
    """
    # Need to compress Rfam.cm for extract_rfam_models.py
    gzip -c ${rfam_cm[0]} > Rfam.cm.gz

    extract_rfam_models.py Rfam.cm.gz ${family_txt_gz} ${rna_types_file} subset.cm selected_families.txt
    """
}

/*
 * Process 3D: Filter clanin file for selected families
 */
process FILTER_CLANIN {
    tag "filter_clanin"
    cpus 1
    memory '1 GB'

    input:
    path clanin
    path families_file

    output:
    path "subset.clanin"

    script:
    """
    filter_rfam_clanin.py ${clanin} ${families_file} subset.clanin
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
    publishDir "${params.outdir}/genomes", mode: 'copy', pattern: "*.fa", optional: true

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
    publishDir "${params.outdir}/cmscan", mode: 'copy', pattern: "*.tblout"
    errorStrategy 'ignore'

    input:
    tuple val(accession), path(genome_fna), path(rfam_cm), path(clanin)

    output:
    tuple val(accession), path("${accession}.tblout"), path(genome_fna)

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
 * Process 6: Extract hit sequences from genome using esl-sfetch
 */
process EXTRACT_HIT_SEQUENCES {
    tag "${accession}"
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/sequences", mode: 'copy', pattern: "*.fa"

    input:
    tuple val(accession), path(tblout), path(genome_fna)

    output:
    tuple val(accession), path("${accession}_hits.fa"), path(tblout)

    script:
    """
    # Parse tblout file and extract sequences using esl-sfetch
    # tblout format (--fmt 2): columns are space-separated
    # Column 1: target name (Rfam accession)
    # Column 2: target accession
    # Column 3: query name (sequence name)
    # Column 4: query accession
    # Columns 8-9: seq from/to (start-stop positions)

    # Create index for genome file
    esl-sfetch --index ${genome_fna}

    # Extract coordinates and Rfam accessions from tblout
    grep -v '^#' ${tblout} | awk '{print \$3"/"\$8"-"\$9}' > coords.txt
    grep -v '^#' ${tblout} | awk '{print \$1}' > rfam_accs.txt

    # Fetch all sequences at once using esl-sfetch
    # -C: use coordinate format (seqname/start-stop)
    # -f: read coordinates from file
    if [ -s coords.txt ]; then
        esl-sfetch -Cf ${genome_fna} coords.txt | \\
        paste -d' ' - rfam_accs.txt | \\
        awk '/^>/ {print \$1, \$2; next} {print}' > ${accession}_hits.fa
    else
        # Create empty file if no hits
        touch ${accession}_hits.fa
    fi
    """
}

/*
 * Process 6B: Identify models from selected families that have no above-threshold hits
 */
process IDENTIFY_MISSING_MODELS {
    tag "${accession}"
    cpus 1
    memory '1 GB'

    input:
    tuple val(accession), path(tblout), path(selected_families)

    output:
    tuple val(accession), path("${accession}.missing_models.txt")

    script:
    """
    identify_missing_models.py ${tblout} ${selected_families} ${accession}.missing_models.txt --genome-id ${accession}
    """
}

/*
 * Process 6C: Extract missing models from full Rfam.cm using cmfetch
 */
process EXTRACT_MISSING_MODELS {
    tag "${accession}"
    cpus 1
    memory '4 GB'

    input:
    tuple val(accession), path(missing_models_file), path(rfam_cm)

    output:
    tuple val(accession), path("${accession}.missing.cm*"), optional: true

    script:
    """
    # Check if there are any missing models
    if [ -s ${missing_models_file} ]; then
        # Extract each missing model using cmfetch
        while read -r model_acc; do
            cmfetch ${rfam_cm[0]} "\${model_acc}" >> ${accession}.missing.cm
        done < ${missing_models_file}

        # Press the CM file for faster searching
        cmpress ${accession}.missing.cm
    else
        echo "No missing models for ${accession}"
    fi
    """
}

/*
 * Process 6D: Run cmscan without threshold on missing models only
 */
process CMSCAN_BELOW_THRESHOLD {
    tag "${accession}"
    cpus 4
    memory '8 GB'
    time '4 h'

    input:
    tuple val(accession), path(missing_cm), path(genome_fna), path(clanin)

    output:
    tuple val(accession), path("${accession}.below_threshold.tblout")

    script:
    """
    # Run cmscan with:
    # -T 0: report all hits with score > 0 (no threshold)
    # --rfam: use Rfam-specific settings
    # --nohmmonly: don't use HMM-only mode for any models
    # --clanin: enable clan competition to avoid overlapping hits
    # --fmt 2: use standard tabular output format
    cmscan -T 20 \\
        --rfam \\
        --nohmmonly \\
        --tblout ${accession}.below_threshold.tblout \\
        --fmt 2 \\
        --clanin ${clanin} \\
        --cpu ${task.cpus} \\
        ${missing_cm[0]} \\
        ${genome_fna}
    """
}

/*
 * Process 6E: Extract best scoring hit per model from below-threshold results
 */
process EXTRACT_BEST_HITS {
    tag "${accession}"
    cpus 1
    memory '2 GB'

    input:
    tuple val(accession), path(tblout)

    output:
    tuple val(accession), path("${accession}.best_hits.tblout")

    script:
    """
    # Use extract_best_hits.py to get the best hit per model
    # Then format as tblout (strip the header comments)
    extract_best_hits.py ${tblout} temp_best_hits.txt --genome-id ${accession}

    # Extract just the hit lines (skip comments) to create a valid tblout file
    grep -v '^#' temp_best_hits.txt > ${accession}.best_hits.tblout || touch ${accession}.best_hits.tblout
    """
}

/*
 * Process 6F: Extract sequences from below-threshold hits
 */
process EXTRACT_HIT_SEQUENCES_BELOW_THRESHOLD {
    tag "${accession}"
    cpus 1
    memory '2 GB'
    publishDir "${params.outdir}/sequences_below_threshold", mode: 'copy', pattern: "*.fa"

    input:
    tuple val(accession), path(tblout), path(genome_fna)

    output:
    tuple val(accession), path("${accession}_below_threshold_hits.fa"), path(tblout)

    script:
    """
    # Parse tblout file and extract sequences using esl-sfetch
    # tblout format (--fmt 2): columns are space-separated
    # Column 1: target name (Rfam accession)
    # Column 2: target accession
    # Column 3: query name (sequence name)
    # Column 4: query accession
    # Columns 8-9: seq from/to (start-stop positions)

    # Create index for genome file
    esl-sfetch --index ${genome_fna}

    # Extract coordinates and Rfam accessions from tblout
    grep -v '^#' ${tblout} | awk '{print \$3"/"\$8"-"\$9}' > coords.txt
    grep -v '^#' ${tblout} | awk '{print \$1}' > rfam_accs.txt

    # Fetch all sequences at once using esl-sfetch
    # -C: use coordinate format (seqname/start-stop)
    # -f: read coordinates from file
    if [ -s coords.txt ]; then
        esl-sfetch -Cf ${genome_fna} coords.txt | \\
        paste -d' ' - rfam_accs.txt | \\
        awk '/^>/ {print \$1, \$2; next} {print}' > ${accession}_below_threshold_hits.fa
    else
        # Create empty file if no hits
        touch ${accession}_below_threshold_hits.fa
    fi
    """
}

/*
 * Process 7: Build SQLite database from all hits and Rfam metadata
 */
process BUILD_SQLITE_DATABASE {
    tag "build_database"
    cpus 2
    memory '8 GB'
    time '2 h'
    publishDir "${params.outdir}/database", mode: 'copy'

    input:
    path above_threshold_tblout_files
    path above_threshold_fasta_files
    path below_threshold_tblout_files
    path below_threshold_fasta_files
    path bacdive_jsonl
    tuple path(family_txt), path(clan_txt), path(clan_membership_txt)

    output:
    path "rfam_hits.db"

    script:
    """
    # Create mapping file: sequence_accession -> genome_id
    echo "Creating sequence to genome mapping..."
    > sequence_to_genome.map
    for tblout in ${above_threshold_tblout_files}; do
        # Extract genome ID from filename (e.g., GCA_000001234.tblout -> GCA_000001234)
        genome_id=\$(basename "\${tblout}" .tblout)
        # Extract unique sequence accessions (column 3) and map to genome_id
        grep -v '^#' "\${tblout}" | awk -v gid="\${genome_id}" '{print \$3 "\\t" gid}' >> sequence_to_genome.map || true
    done
    for tblout in ${below_threshold_tblout_files}; do
        genome_id=\$(basename "\${tblout}" .best_hits.tblout)
        genome_id=\$(basename "\${genome_id}" .tblout)
        grep -v '^#' "\${tblout}" | awk -v gid="\${genome_id}" '{print \$3 "\\t" gid}' >> sequence_to_genome.map || true
    done
    # Remove duplicates and sort
    sort -u sequence_to_genome.map > sequence_to_genome_unique.map

    # Concatenate all above-threshold tblout files (strip comments)
    echo "Concatenating above-threshold tblout files..."
    for tblout in ${above_threshold_tblout_files}; do
        grep -v '^#' "\${tblout}" >> above_threshold_all.tblout || true
    done

    # Concatenate all above-threshold FASTA files
    echo "Concatenating above-threshold FASTA files..."
    cat ${above_threshold_fasta_files} > above_threshold_all.fa

    # Concatenate all below-threshold tblout files (strip comments)
    echo "Concatenating below-threshold tblout files..."
    for tblout in ${below_threshold_tblout_files}; do
        grep -v '^#' "\${tblout}" >> below_threshold_all.tblout || true
    done

    # Concatenate all below-threshold FASTA files
    echo "Concatenating below-threshold FASTA files..."
    cat ${below_threshold_fasta_files} > below_threshold_all.fa

    # Build database with concatenated files
    build_sqlite_db.py rfam_hits.db \\
        --above-tblout above_threshold_all.tblout \\
        --above-fasta above_threshold_all.fa \\
        --below-tblout below_threshold_all.tblout \\
        --below-fasta below_threshold_all.fa \\
        --sequence-map sequence_to_genome_unique.map \\
        --bacdive-jsonl ${bacdive_jsonl} \\
        --family-txt ${family_txt} \\
        --clan-txt ${clan_txt} \\
        --clan-membership-txt ${clan_membership_txt}
    """
}

/*
 * Process 8: Create alignments for all families with hits
 */
process CREATE_ALIGNMENTS {
    tag "create_alignments"
    cpus 4
    memory '16 GB'
    time '12 h'
    publishDir "${params.outdir}/database", mode: 'copy', pattern: "*.db"
    publishDir "${params.outdir}/alignments", mode: 'copy', pattern: "*.sto"

    input:
    path database
    path rfam_cm

    output:
    path "rfam_hits_with_alignments.db"
    path "*.sto", optional: true

    script:
    """
    # Copy database to output name
    cp ${database} rfam_hits_with_alignments.db

    # Create alignments and store in database, saving Stockholm files to current directory
    create_alignments.py rfam_hits_with_alignments.db ${rfam_cm[0]} --output-dir .
    """
}

/*
 * Workflow
 */
workflow {
    // Convert CSV to JSONL
    Channel.fromPath(params.input_csv, checkIfExists: true) | set { bacdrive }

    // Download Rfam covariance models
    DOWNLOAD_RFAM | set { rfam }

    // Download Rfam metadata for database building
    DOWNLOAD_RFAM_METADATA | set { rfam_metadata }

    // Load RNA types file for below-threshold scanning
    Channel.fromPath(params.rna_types_file, checkIfExists: true) | set { rna_types }

    // Extract selected RNA type models and filter clanin
    EXTRACT_MODELS_BY_TYPE(rfam, rfam_metadata, rna_types) | set { models_and_families }

    // Split the tuple to get files separately
    models_and_families \
    | map { cm_files, families -> families } \
    | set { selected_families }

    // Get the first element from rfam tuple for filtering clanin
    rfam \
    | map { cm_files, clanin -> clanin } \
    | set { rfam_clanin }

    FILTER_CLANIN(rfam_clanin, selected_families) | set { subset_clanin }

    // Process genomes
    bacdrive \
    | CSV_TO_JSONL \
    | set { bacdive_jsonl }

    bacdive_jsonl \
    | EXTRACT_GCA_ACCESSIONS \
    | splitText \
    | map { it.trim() } \
    | filter { it.length() > 0 } \
    | DOWNLOAD_GENOME \
    | set { genomes }

    // Run above-threshold cmscan and extract sequences
    genomes \
    | combine(rfam) \
    | CMSCAN \
    | EXTRACT_HIT_SEQUENCES \
    | set { above_threshold_results }

    // Identify missing models for each genome
    above_threshold_results \
    | map { accession, fasta, tblout -> [accession, tblout] } \
    | combine(selected_families) \
    | IDENTIFY_MISSING_MODELS \
    | set { missing_models_list }

    // Extract missing models from full Rfam.cm
    missing_models_list \
    | combine(rfam.map { cm_files, clanin -> cm_files }) \
    | EXTRACT_MISSING_MODELS \
    | set { missing_cm_files }

    // Combine missing CM files with genome and clanin for below-threshold scan
    // Need to join back with genomes by accession
    missing_cm_files \
    | join(genomes) \
    | map { accession, missing_cm, genome_fna -> [accession, missing_cm, genome_fna] } \
    | combine(subset_clanin) \
    | CMSCAN_BELOW_THRESHOLD \
    | EXTRACT_BEST_HITS \
    | set { best_hits_results }

    // Join best hits results with genomes to get genome FASTA for sequence extraction
    best_hits_results \
    | join(genomes) \
    | map { accession, best_tblout, genome_fna -> [accession, best_tblout, genome_fna] } \
    | EXTRACT_HIT_SEQUENCES_BELOW_THRESHOLD \
    | set { below_threshold_results }

    // Collect above-threshold results
    above_threshold_results \
    | map { accession, fasta, tblout -> tblout } \
    | collect \
    | set { above_tblouts }

    above_threshold_results \
    | map { accession, fasta, tblout -> fasta } \
    | collect \
    | set { above_fastas }

    // Collect below-threshold results
    below_threshold_results \
    | map { accession, fasta, tblout -> tblout } \
    | collect \
    | set { below_tblouts }

    below_threshold_results \
    | map { accession, fasta, tblout -> fasta } \
    | collect \
    | set { below_fastas }

    // Build SQLite database with both above and below threshold hits
    BUILD_SQLITE_DATABASE(above_tblouts, above_fastas, below_tblouts, below_fastas, bacdive_jsonl, rfam_metadata) \
    | set { database }

    // Create alignments for all families and add to database
    CREATE_ALIGNMENTS(database, rfam.map { cm_files, clanin -> cm_files })
}
