#!/usr/bin/env python3
"""
Build SQLite database from cmscan results and Rfam metadata.

This script creates a comprehensive SQLite database containing:
- Hit information from cmscan tblout files
- Hit sequences extracted from genomes
- Rfam family metadata
- Rfam clan metadata
- Clan membership information
- BacDive genome metadata (temperature, etc.)
"""

import json
from pathlib import Path
from typing import Dict, List

import click
from sqlite_utils import Database


def parse_tblout_line(line: str) -> Dict[str, any]:
    """Parse a single line from cmscan tblout output (format 2).

    Returns a dictionary with all tblout columns.
    Format 2 includes an idx column as the first column.
    """
    parts = line.strip().split()

    # Infernal cmscan tblout format 2 has idx + 18 standard columns
    # (plus optional description which may contain spaces)
    if len(parts) < 19:
        return None

    return {
        "idx": int(parts[0]),
        "target_name": parts[1],  # Family ID (e.g., "5S_rRNA")
        "rfam_acc": parts[2],  # Family accession (e.g., "RF00001") - FK to family.rfam_acc
        "query_name": parts[3],
        "query_accession": parts[4] if parts[4] != "-" else None,
        "clan_name": parts[5] if parts[5] != "-" else None,
        "mdl": parts[6],
        "mdl_from": int(parts[7]),
        "mdl_to": int(parts[8]),
        "seq_from": int(parts[9]),
        "seq_to": int(parts[10]),
        "strand": parts[11],
        "trunc": parts[12],
        "pass_": parts[13],  # 'pass' is a Python keyword, use 'pass_'
        "gc": float(parts[14]),
        "bias": float(parts[15]),
        "score": float(parts[16]),
        "e_value": float(parts[17]) if parts[17] != "-" else None,
        "inc": parts[18],
        "olp": parts[19] if len(parts) > 19 and parts[19] != "-" else None,
        "anyidx": int(parts[20]) if len(parts) > 20 and parts[20] != "-" else None,
        "afrct1": float(parts[21]) if len(parts) > 21 and parts[21] != "-" and parts[21] != "\"" else None,
        "afrct2": float(parts[22]) if len(parts) > 22 and parts[22] != "-" and parts[22] != "\"" else None,
        "winidx": parts[23] if len(parts) > 23 and parts[23] != "-" and parts[23] != "\"" else None,
        "wfrct1": float(parts[24]) if len(parts) > 24 and parts[24] != "-" and parts[24] != "\"" else None,
        "wfrct2": float(parts[25]) if len(parts) > 25 and parts[25] != "-" and parts[25] != "\"" else None,
        "mdl_len": int(parts[26]) if len(parts) > 26 else None,
        "seq_len": int(parts[27]) if len(parts) > 27 else None,
        "description": " ".join(parts[28:]) if len(parts) > 28 else None,
    }


def parse_fasta(fasta_path: Path) -> Dict[str, str]:
    """Parse FASTA file and return dict of {header: sequence}.

    Headers are stored without the leading '>'.
    """
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # Save previous sequence
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Save last sequence
        if current_header:
            sequences[current_header] = "".join(current_seq)

    return sequences


# Schema definitions from Rfam SQL files
# https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/
RFAM_SCHEMAS = {
    'family': [
        'rfam_acc', 'rfam_id', 'auto_wiki', 'description', 'author', 'seed_source',
        'gathering_cutoff', 'trusted_cutoff', 'noise_cutoff', 'comment', 'previous_id',
        'cmbuild', 'cmcalibrate', 'cmsearch', 'num_seed', 'num_full', 'num_genome_seq',
        'num_refseq', 'type', 'structure_source', 'number_of_species', 'number_3d_structures',
        'num_pseudonokts', 'tax_seed', 'ecmli_lambda', 'ecmli_mu', 'ecmli_cal_db',
        'ecmli_cal_hits', 'maxl', 'clen', 'match_pair_node', 'hmm_tau', 'hmm_lambda',
        'created', 'updated'
    ],
    'clan': [
        'clan_acc', 'id', 'previous_id', 'description', 'author', 'comment', 'created', 'updated'
    ],
    'clan_membership': [
        'clan_acc', 'rfam_acc'
    ]
}


def parse_rfam_table(txt_path: Path, table_name: str) -> List[Dict]:
    """Parse tab-delimited Rfam table file using proper schema.

    Args:
        txt_path: Path to the tab-delimited text file
        table_name: Name of the table (must match RFAM_SCHEMAS key)

    Returns:
        List of dictionaries with column names as keys.
        MySQL dumps use \\N to represent NULL values which are converted to None.
    """
    if table_name not in RFAM_SCHEMAS:
        raise ValueError(f"Unknown table '{table_name}'. Must be one of: {list(RFAM_SCHEMAS.keys())}")

    header = RFAM_SCHEMAS[table_name]
    rows = []

    with open(txt_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')

            # Create row dict matching schema
            row = {}
            for i, col in enumerate(header):
                if i < len(parts):
                    # Convert MySQL NULL marker \N to Python None
                    value = parts[i]
                    row[col] = None if value == '\\N' else value
                else:
                    # Missing columns default to None
                    row[col] = None
            rows.append(row)

    return rows


@click.command()
@click.argument('output_db', type=click.Path())
@click.option('--above-tblout', type=click.Path(exists=True), help='Above-threshold tblout file (concatenated)')
@click.option('--above-fasta', type=click.Path(exists=True), help='Above-threshold FASTA file (concatenated)')
@click.option('--below-tblout', type=click.Path(exists=True), help='Below-threshold tblout file (concatenated)')
@click.option('--below-fasta', type=click.Path(exists=True), help='Below-threshold FASTA file (concatenated)')
@click.option('--sequence-map', required=True, type=click.Path(exists=True), help='Sequence accession to genome ID mapping file')
@click.option('--bacdive-jsonl', required=True, type=click.Path(exists=True), help='BacDive JSONL file with genome metadata')
@click.option('--family-txt', required=True, type=click.Path(exists=True), help='family.txt file')
@click.option('--clan-txt', required=True, type=click.Path(exists=True), help='clan.txt file')
@click.option('--clan-membership-txt', required=True, type=click.Path(exists=True), help='clan_membership.txt file')
def main(output_db, above_tblout, above_fasta, below_tblout, below_fasta, sequence_map, bacdive_jsonl, family_txt, clan_txt, clan_membership_txt):
    """Build SQLite database from cmscan results (above and below threshold) and Rfam metadata."""

    click.echo(f"Creating database: {output_db}")
    db = Database(output_db)

    # Enable foreign key enforcement in SQLite
    db.execute("PRAGMA foreign_keys = ON")
    click.echo("Foreign key enforcement enabled")

    # Load sequence accession to genome ID mapping
    click.echo("Loading sequence to genome mapping...")
    seq_to_genome = {}
    with open(sequence_map) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                seq_acc, genome_id = parts
                seq_to_genome[seq_acc] = genome_id
    click.echo(f"  Loaded {len(seq_to_genome)} sequence accessions")

    # Create and populate sequences table
    click.echo("Creating sequences table...")
    sequences_data = [
        {"accession": acc, "genome_id": genome_id}
        for acc, genome_id in seq_to_genome.items()
    ]
    db["sequences"].insert_all(sequences_data, pk="accession", replace=True, strict=True)
    click.echo(f"  Inserted {len(sequences_data)} sequences")

    # Load BacDive data and split into genomes, strains, and measurements
    click.echo("Loading BacDive data...")
    genomes_dict = {}  # genome_id -> genome record
    strains = []  # List of strain records
    measurements = []  # List of measurement records
    strain_id_counter = 1

    with open(bacdive_jsonl) as f:
        for line in f:
            record = json.loads(line)

            # Extract genome information
            genome_data = record.get("genome", {})
            genome_acc = genome_data.get("accession", "").strip()

            # Only process records with valid GCA accessions
            if genome_acc.startswith("GCA_"):
                # Add genome to dict (will deduplicate automatically)
                if genome_acc not in genomes_dict:
                    genomes_dict[genome_acc] = {
                        "genome_id": genome_acc,
                        "genome_sequence_database": genome_data.get("database", "")
                    }

                # Extract strain information
                strain_data = record.get("strain", {})
                strain_record = {
                    "strain_id": strain_id_counter,
                    "genome_id": genome_acc,
                    "bacdive_id": strain_data.get("bacdive_id", ""),
                    "designation": strain_data.get("designation", ""),
                    "strain_numbers": strain_data.get("strain_numbers", ""),
                    "is_type_strain": strain_data.get("is_type_strain", False),
                    "species": strain_data.get("species", "")
                }
                strains.append(strain_record)

                # Extract measurement information
                measurement_data = record.get("measurement")
                if measurement_data:
                    measurement_record = {
                        "strain_id": strain_id_counter,
                        "property_type": measurement_data.get("property_type", ""),
                        "measurement_kind": measurement_data.get("measurement_kind", ""),
                        "value_min": measurement_data.get("value_min"),
                        "value_max": measurement_data.get("value_max"),
                        "value_unit": measurement_data.get("value_unit", ""),
                        "test_result": measurement_data.get("test_result", "")
                    }
                    measurements.append(measurement_record)

                strain_id_counter += 1

    # Insert genomes
    genomes_list = list(genomes_dict.values())
    db["genomes"].insert_all(genomes_list, pk="genome_id", replace=True, strict=True)
    click.echo(f"  Inserted {len(genomes_list)} unique genomes")

    # Insert strains
    db["strains"].insert_all(strains, pk="strain_id", replace=True, strict=True)
    click.echo(f"  Inserted {len(strains)} strains")

    # Insert measurements
    if measurements:
        db["measurements"].insert_all(measurements, replace=True, strict=True)
        click.echo(f"  Inserted {len(measurements)} measurements")

    # Parse Rfam metadata files
    click.echo("Loading Rfam family metadata...")
    families = parse_rfam_table(Path(family_txt), 'family')
    db["family"].insert_all(families, pk="rfam_acc", replace=True, strict=True)
    click.echo(f"  Inserted {len(families)} families")

    click.echo("Loading Rfam clan metadata...")
    clans = parse_rfam_table(Path(clan_txt), 'clan')
    db["clan"].insert_all(clans, pk="clan_acc", replace=True, strict=True)
    click.echo(f"  Inserted {len(clans)} clans")

    click.echo("Loading Rfam clan membership...")
    memberships = parse_rfam_table(Path(clan_membership_txt), 'clan_membership')
    db["clan_membership"].insert_all(memberships, replace=True, strict=True)
    click.echo(f"  Inserted {len(memberships)} clan memberships")

    # Add foreign key constraints for clan_membership
    click.echo("Adding foreign key constraints for clan_membership...")
    db["clan_membership"].add_foreign_key("clan_acc", "clan", "clan_acc", ignore=True)
    db["clan_membership"].add_foreign_key("rfam_acc", "family", "rfam_acc", ignore=True)
    click.echo("  Foreign keys added")

    # Helper function to process a single FASTA and tblout file
    def process_hits(tblout_file: Path, fasta_file: Path, hit_type: str):
        """Process hits from concatenated tblout and fasta files, returning count of hits inserted."""
        # Parse FASTA file and build sequence lookup
        click.echo(f"Loading {hit_type} hit sequences from FASTA file...")
        all_sequences = {}

        sequences = parse_fasta(fasta_file)
        click.echo(f"  Loaded {len(sequences)} sequences")

        # Store sequences with lookup key
        for header, seq in sequences.items():
            # Header format: "seqname/start-stop"
            # Example: "CP001878.2/12345-12456"
            # Use the entire header as the key (no additional fields)
            coord_part = header.split()[0] if ' ' in header else header
            all_sequences[coord_part] = seq

        # Parse tblout file and create hit records
        click.echo(f"Loading {hit_type} hit data from tblout file...")
        hit_count = 0
        genomes_processed = set()
        missing_mapping_count = 0

        with open(tblout_file) as f:
            for line in f:
                # Skip comment lines (should already be stripped, but be safe)
                if line.startswith('#'):
                    continue

                # Skip empty lines
                if not line.strip():
                    continue

                hit = parse_tblout_line(line)
                if not hit:
                    continue

                # Look up genome_id from sequence accession (query_name)
                query_name = hit['query_name']
                genome_id = seq_to_genome.get(query_name)

                if not genome_id:
                    missing_mapping_count += 1
                    if missing_mapping_count <= 5:  # Only show first 5 warnings
                        click.echo(f"  Warning: No genome mapping for sequence {query_name}", err=True)
                    continue

                genomes_processed.add(genome_id)

                # Add genome_id and hit_type to the hit record
                hit['genome_id'] = genome_id
                hit['hit_type'] = hit_type

                # Look up the sequence for this hit
                # Create coordinate key: query_name/seq_from-seq_to
                coord_key = f"{query_name}/{hit['seq_from']}-{hit['seq_to']}"

                hit['hit_sequence'] = all_sequences.get(coord_key, None)

                if hit['hit_sequence'] is None:
                    click.echo(f"  Warning: No sequence found for {coord_key}", err=True)

                # Insert hit into database
                db["hits"].insert(hit, strict=True)
                hit_count += 1

        if missing_mapping_count > 5:
            click.echo(f"  Warning: {missing_mapping_count} total sequences without genome mapping", err=True)

        click.echo(f"  Processed {len(genomes_processed)} genomes")
        return hit_count

    # Process above-threshold hits
    above_count = 0
    if above_tblout and above_fasta:
        above_count = process_hits(Path(above_tblout), Path(above_fasta), "above_threshold")
        click.echo(f"Total above-threshold hits inserted: {above_count}")

    # Process below-threshold hits
    below_count = 0
    if below_tblout and below_fasta:
        below_count = process_hits(Path(below_tblout), Path(below_fasta), "below_threshold")
        click.echo(f"Total below-threshold hits inserted: {below_count}")

    click.echo(f"Total hits inserted: {above_count + below_count}")

    # Create indexes for better query performance
    click.echo("Creating indexes...")
    db["hits"].create_index(["genome_id"], if_not_exists=True)
    db["hits"].create_index(["rfam_acc"], if_not_exists=True)
    db["hits"].create_index(["query_name"], if_not_exists=True)
    db["hits"].create_index(["e_value"], if_not_exists=True)
    db["hits"].create_index(["score"], if_not_exists=True)
    db["hits"].create_index(["hit_type"], if_not_exists=True)

    # Create foreign key relationships (informational, not enforced by default)
    click.echo("Adding foreign key references...")
    try:
        # sequences.genome_id -> genomes.genome_id
        db["sequences"].add_foreign_key("genome_id", "genomes", "genome_id")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for sequences->genomes: {e}", err=True)

    try:
        # hits.rfam_acc -> family.rfam_acc
        db["hits"].add_foreign_key("rfam_acc", "family", "rfam_acc")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->family: {e}", err=True)

    try:
        # hits.query_name -> sequences.accession
        db["hits"].add_foreign_key("query_name", "sequences", "accession")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->sequences: {e}", err=True)

    try:
        # hits.genome_id -> genomes.genome_id
        db["hits"].add_foreign_key("genome_id", "genomes", "genome_id")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->genomes: {e}", err=True)

    try:
        # hits.clan_name -> clan.clan_acc
        db["hits"].add_foreign_key("clan_name", "clan", "clan_acc")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->clan: {e}", err=True)

    try:
        # hits.target_name -> family.rfam_id
        db["hits"].add_foreign_key("target_name", "family", "rfam_id")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->family(rfam_id): {e}", err=True)

    try:
        # strains.genome_id -> genomes.genome_id
        db["strains"].add_foreign_key("genome_id", "genomes", "genome_id")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for strains->genomes: {e}", err=True)

    try:
        # measurements.strain_id -> strains.strain_id
        db["measurements"].add_foreign_key("strain_id", "strains", "strain_id")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for measurements->strains: {e}", err=True)

    # Create indexes for strains and measurements
    click.echo("Creating indexes for strains and measurements...")
    db["strains"].create_index(["genome_id"], if_not_exists=True)
    db["strains"].create_index(["species"], if_not_exists=True)
    db["measurements"].create_index(["strain_id"], if_not_exists=True)
    db["measurements"].create_index(["property_type"], if_not_exists=True)
    db["measurements"].create_index(["measurement_kind"], if_not_exists=True)
    db["measurements"].create_index(["value_min"], if_not_exists=True)
    db["measurements"].create_index(["value_max"], if_not_exists=True)

    # Create empty alignment tables (will be populated by create_alignments.py)
    click.echo("Creating alignment table schemas...")

    # Create alignments table
    db["alignments"].create({
        "alignment_id": int,
        "rfam_acc": str,  # Renamed from family_id for consistency
        "genome_id": str,
        "hit_id": str,  # Reference to hit (accession/seq_from-seq_to)
        "aligned_sequence": str,
        "alignment_start": int,
        "alignment_end": int,
    }, pk="alignment_id", if_not_exists=True, strict=True)

    # Create alignment_metadata table
    db["alignment_metadata"].create({
        "rfam_acc": str,  # Renamed from family_id for consistency
        "num_sequences": int,
        "alignment_length": int,
        "consensus_structure": str,
        "created_date": str,
    }, pk="rfam_acc", if_not_exists=True, strict=True)

    # Create alignment_files table to store complete Stockholm alignments as blobs
    db["alignment_files"].create({
        "rfam_acc": str,  # Renamed from family_id for consistency
        "stockholm_alignment": bytes,  # BLOB of complete Stockholm file
        "file_size": int,
        "created_date": str,
    }, pk="rfam_acc", if_not_exists=True, strict=True)

    click.echo("  Created alignments, alignment_metadata, and alignment_files tables (empty)")

    # Create SQL view for best hits (de-overlapped)
    click.echo("Creating best_hits view (excluding overlapping hits)...")
    db.execute("""
        CREATE VIEW IF NOT EXISTS best_hits AS
        SELECT *
        FROM hits
        WHERE olp NOT LIKE '%=%'
    """)
    click.echo("  Created best_hits view (filters out lower-scoring overlapping hits)")

    # Create SQL view for aggregated genome properties
    click.echo("Creating genome_properties_view for aggregated data...")
    db.execute("""
        CREATE VIEW IF NOT EXISTS genome_properties_view AS
        SELECT
            g.genome_id,
            s.species,
            MIN(CASE WHEN m.property_type = 'temperature' AND m.measurement_kind = 'growth' THEN m.value_min END) as growth_temp_min,
            MAX(CASE WHEN m.property_type = 'temperature' AND m.measurement_kind = 'growth' THEN m.value_max END) as growth_temp_max,
            AVG(CASE WHEN m.property_type = 'temperature' AND m.measurement_kind = 'optimum' THEN m.value_min END) as optimal_temp,
            COUNT(CASE WHEN m.property_type = 'temperature' THEN 1 END) as temp_measurement_count
        FROM genomes g
        LEFT JOIN strains s ON g.genome_id = s.genome_id
        LEFT JOIN measurements m ON s.strain_id = m.strain_id
        GROUP BY g.genome_id, s.species
    """)
    click.echo("  Created genome_properties_view")

    click.echo(f"\nDatabase created successfully: {output_db}")
    click.echo(f"  Tables: {', '.join(db.table_names())}")

    # Print example queries
    click.echo("\nExample queries:")
    click.echo("  -- Use best_hits view for genome annotation (excludes overlapping hits):")
    click.echo("  SELECT * FROM best_hits WHERE genome_id = 'GCA_000005825';")
    click.echo("\n  -- Find riboswitch hits in genomes that grow below 20°C (using best hits):")
    click.echo("  SELECT h.*, gp.* FROM best_hits h")
    click.echo("  JOIN genome_properties_view gp ON h.genome_id = gp.genome_id")
    click.echo("  WHERE h.target_name LIKE '%riboswitch%' AND gp.growth_temp_max <= 20;")
    click.echo("\n  -- Count best hits by genome with temperature data:")
    click.echo("  SELECT gp.genome_id, gp.species, gp.growth_temp_min, gp.growth_temp_max,")
    click.echo("         COUNT(*) as hit_count FROM best_hits h")
    click.echo("  JOIN genome_properties_view gp ON h.genome_id = gp.genome_id")
    click.echo("  GROUP BY gp.genome_id, gp.species, gp.growth_temp_min, gp.growth_temp_max;")
    click.echo("\n  -- Compare all hits vs best hits counts:")
    click.echo("  SELECT")
    click.echo("    (SELECT COUNT(*) FROM hits) as total_hits,")
    click.echo("    (SELECT COUNT(*) FROM best_hits) as best_hits,")
    click.echo("    (SELECT COUNT(*) FROM hits WHERE olp LIKE '%=%') as overlapping_hits;")
    click.echo("\n  -- View all temperature measurements for a genome:")
    click.echo("  SELECT s.species, m.* FROM measurements m")
    click.echo("  JOIN strains s ON m.strain_id = s.strain_id")
    click.echo("  WHERE s.genome_id = 'GCA_000005825';")
    click.echo("\n  -- Check for inconsistencies (genomes with >3 temperature measurements):")
    click.echo("  SELECT * FROM genome_properties_view WHERE temp_measurement_count > 3;")
    click.echo("\n  -- View alignment for a specific family:")
    click.echo("  SELECT * FROM alignments WHERE rfam_acc = 'RF00001';")


if __name__ == '__main__':
    main()
