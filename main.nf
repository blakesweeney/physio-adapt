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
    tuple path("family.txt"), path("clan.txt"), path("clan_membership.txt")

    script:
    """
    # Download and uncompress Rfam database files from FTP
    curl -L https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/database_files/family.txt.gz | gunzip > family.txt
    curl -L https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/database_files/clan.txt.gz | gunzip > clan.txt
    curl -L https://ftp.ebi.ac.uk/pub/databases/Rfam/15.0/database_files/clan_membership.txt.gz | gunzip > clan_membership.txt
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
    tuple path(rfam_cm), path(clanin), path(family_txt), path(clan_txt), path(clan_membership_txt), path(rna_types_file)

    output:
    tuple path("subset.cm*"), path("selected_families.txt")

    script:
    """
    # Extract Rfam accessions for families matching RNA types
    # Column 1 is rfam_acc, column 19 is type (e.g., "Cis-reg; riboswitch;")

    # Read RNA types from file (skip comments and empty lines)
    grep -v '^#' ${rna_types_file} | grep -v '^\$' > rna_types_clean.txt

    # Extract column 1 (rfam_acc) and column 19 (type) from family.txt
    # Then grep for lines matching any RNA type, then extract just the accession
    cat ${family_txt} | awk 'BEGIN{FS="\\t"; OFS="\\t"} NR>1 {print \$1, \$19}' | \\
        grep -f rna_types_clean.txt | \\
        awk '{print \$1}' | \\
        sort -u > selected_families.txt

    # Check if we found any families
    NUM_FAMILIES=\$(wc -l < selected_families.txt)

    if [ "\${NUM_FAMILIES}" -eq 0 ]; then
        echo "ERROR: No families found for specified RNA types" >&2
        exit 1
    fi

    # Extract models using cmfetch (it can read a file of accessions with -f)
    cmfetch -f ${rfam_cm[0]} selected_families.txt > subset.cm

    # Press the CM file for faster searching
    cmpress -F subset.cm
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
    tuple path(clanin), path(families_file)

    output:
    path "subset.clanin"

    script:
    """
    # Filter clanin to keep only clans containing at least one selected family
    # Clanin format: CL00001 RF00001 RF00002 RF00003 ...
    # Keep entire line if any family accession matches our selected families

    grep -f ${families_file} ${clanin} > subset.clanin || touch subset.clanin
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
    # Column 3: query name (sequence accession, e.g., GG693512.1)
    # Columns 10-11: seq from/to (start-stop positions)

    # Create index for genome file
    esl-sfetch --index ${genome_fna}

    # Extract coordinates from tblout: query_name/seq_from-seq_to
    grep -v '^#' ${tblout} | awk '{print \$3"/"\$10"-"\$11}' > coords.txt

    # Fetch all sequences at once using esl-sfetch
    # -C: use coordinate format (seqname/start-stop)
    # -f: read coordinates from file
    if [ -s coords.txt ]; then
        esl-sfetch -Cf ${genome_fna} coords.txt > ${accession}_hits.fa
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
    # -T 20: report all hits with score > 0 (no threshold)
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
    # Column 3: query name (sequence accession, e.g., GG693512.1)
    # Columns 10-11: seq from/to (start-stop positions)

    # Create index for genome file
    esl-sfetch --index ${genome_fna}

    # Extract coordinates from tblout: query_name/seq_from-seq_to
    grep -v '^#' ${tblout} | awk '{print \$3"/"\$10"-"\$11}' > coords.txt

    # Fetch all sequences at once using esl-sfetch
    # -C: use coordinate format (seqname/start-stop)
    # -f: read coordinates from file
    if [ -s coords.txt ]; then
        esl-sfetch -Cf ${genome_fna} coords.txt > ${accession}_below_threshold_hits.fa
    else
        # Create empty file if no hits
        touch ${accession}_below_threshold_hits.fa
    fi
    """
}

/*
 * Process 7A: Concatenate above-threshold tblout files and create sequence mapping
 */
process CONCATENATE_ABOVE_THRESHOLD_TBLOUT {
    tag "concat_above_tblout"
    cpus 1
    memory '4 GB'
    time '1 h'

    input:
    path tblout_files

    output:
    tuple path("above_threshold_all.tblout"), path("sequence_to_genome.map")

    script:
    """
    # Concatenate tblout files (strip comments)
    find . -name '*.tblout' -type f -exec grep -v '^#' {} \\; > above_threshold_all.tblout

    # Create sequence to genome mapping
    # For each tblout file, extract filename and column 4 (query_name)
    find . -name '*.tblout' -type f -exec sh -c '
        genome_id=\$(basename "{}" .tblout)
        grep -v "^#" "{}" | awk -v gid="\$genome_id" "{print \\\$4 \"\\t\" gid}"
    ' \\; > sequence_to_genome.map
    """
}

/*
 * Process 7B: Concatenate above-threshold FASTA files
 */
process CONCATENATE_ABOVE_THRESHOLD_FASTA {
    tag "concat_above_fasta"
    cpus 1
    memory '4 GB'
    time '1 h'

    input:
    path fasta_files

    output:
    path "above_threshold_all.fa"

    script:
    """
    find . -name '*.fa' -type f -print0 | xargs -0 cat > above_threshold_all.fa
    """
}

/*
 * Process 7C: Concatenate below-threshold tblout files and add to sequence mapping
 */
process CONCATENATE_BELOW_THRESHOLD_TBLOUT {
    tag "concat_below_tblout"
    cpus 1
    memory '4 GB'
    time '1 h'

    input:
    path tblout_files
    path above_mapping

    output:
    tuple path("below_threshold_all.tblout"), path("sequence_to_genome_all.map")

    script:
    """
    # Concatenate tblout files (strip comments)
    find . -name '*.tblout' -type f -exec grep -v '^#' {} \\; > below_threshold_all.tblout

    # Start with above-threshold mapping and add below-threshold sequences
    cp ${above_mapping} sequence_to_genome_all.map

    # Add below-threshold mappings (handle various filename patterns)
    find . -name '*.tblout' -type f -exec sh -c '
        fname=\$(basename "{}")
        genome_id=\${fname%.below_threshold.tblout}
        genome_id=\${genome_id%.best_hits.tblout}
        genome_id=\${genome_id%.tblout}
        grep -v "^#" "{}" | awk -v gid="\$genome_id" "{print \\\$4 \"\\t\" gid}"
    ' \\; >> sequence_to_genome_all.map

    # Remove duplicates and sort
    sort -u sequence_to_genome_all.map -o sequence_to_genome_all.map
    """
}

/*
 * Process 7D: Concatenate below-threshold FASTA files
 */
process CONCATENATE_BELOW_THRESHOLD_FASTA {
    tag "concat_below_fasta"
    cpus 1
    memory '4 GB'
    time '1 h'

    input:
    path fasta_files

    output:
    path "below_threshold_all.fa"

    script:
    """
    find . -name '*.fa' -type f -print0 | xargs -0 cat > below_threshold_all.fa
    """
}

/*
 * Process 7E: Build SQLite database from concatenated files
 */
process BUILD_SQLITE_DATABASE {
    tag "build_database"
    cpus 2
    memory '8 GB'
    time '2 h'
    publishDir "${params.outdir}/database", mode: 'copy'

    input:
    tuple path(above_tblout), path(seq_map_partial), path(above_fasta), path(below_tblout), path(seq_map_complete), path(below_fasta), path(bacdive_jsonl), path(family_txt), path(clan_txt), path(clan_membership_txt)

    output:
    path "rfam_hits.db"

    script:
    """
    # Build database with concatenated files
    build_sqlite_db.py rfam_hits.db \\
        --above-tblout ${above_tblout} \\
        --above-fasta ${above_fasta} \\
        --below-tblout ${below_tblout} \\
        --below-fasta ${below_fasta} \\
        --sequence-map ${seq_map_complete} \\
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
    tuple path(database), path(rfam_cm)

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
    rfam \
    | combine(rfam_metadata) \
    | combine(rna_types) \
    | EXTRACT_MODELS_BY_TYPE \
    | set { models_and_families }

    // Split the tuple to get files separately
    models_and_families \
    | map { cm_files, families -> families } \
    | set { selected_families }

    // Get the first element from rfam tuple for filtering clanin
    rfam \
    | map { cm_files, clanin -> clanin } \
    | set { rfam_clanin }

    rfam_clanin \
    | combine(selected_families) \
    | FILTER_CLANIN \
    | set { subset_clanin }

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

    // Concatenate above-threshold results
    above_threshold_results \
    | map { accession, fasta, tblout -> tblout } \
    | collect \
    | CONCATENATE_ABOVE_THRESHOLD_TBLOUT \
    | set { above_concat_tblout_and_map }

    above_threshold_results \
    | map { accession, fasta, tblout -> fasta } \
    | collect \
    | CONCATENATE_ABOVE_THRESHOLD_FASTA \
    | set { above_concat_fasta }

    // Concatenate below-threshold results
    below_threshold_results \
    | map { accession, fasta, tblout -> tblout } \
    | collect \
    | combine(above_concat_tblout_and_map.map { tblout, mapping -> mapping }) \
    | CONCATENATE_BELOW_THRESHOLD_TBLOUT \
    | set { below_concat_tblout_and_map }

    below_threshold_results \
    | map { accession, fasta, tblout -> fasta } \
    | collect \
    | CONCATENATE_BELOW_THRESHOLD_FASTA \
    | set { below_concat_fasta }

    // Build SQLite database with concatenated files
    above_concat_tblout_and_map \
    | combine(above_concat_fasta) \
    | combine(below_concat_tblout_and_map) \
    | combine(below_concat_fasta) \
    | combine(bacdive_jsonl) \
    | combine(rfam_metadata) \
    | BUILD_SQLITE_DATABASE \
    | set { database }

    // Create alignments for all families and add to database
    database \
    | combine(rfam.map { cm_files, clanin -> cm_files }) \
    | CREATE_ALIGNMENTS
}
