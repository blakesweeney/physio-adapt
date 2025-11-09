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
        "target_name": parts[1],
        "target_accession": parts[2],
        "query_name": parts[3],
        "query_accession": parts[4],
        "clan_name": parts[5],
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
        "e_value": float(parts[17]),
        "inc": parts[18],
        "olp": parts[19] if len(parts) > 19 else None,
        "anyidx": int(parts[20]) if len(parts) > 20 else None,
        "afrct1": float(parts[21]) if len(parts) > 21 else None,
        "afrct2": float(parts[22]) if len(parts) > 22 else None,
        "winidx": parts[23] if len(parts) > 23 else None,
        "wfrct1": float(parts[24]) if len(parts) > 24 else None,
        "wfrct2": float(parts[25]) if len(parts) > 25 else None,
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


def parse_rfam_table(txt_path: Path) -> List[Dict]:
    """Parse tab-delimited Rfam table file.

    Returns list of dictionaries with column names as keys.
    """
    rows = []

    with open(txt_path, 'r') as f:
        # First line is header
        header = next(f).strip().split('\t')

        for line in f:
            parts = line.strip().split('\t')
            # Handle lines with fewer columns than header
            # (fill missing values with None)
            row = {}
            for i, col in enumerate(header):
                row[col] = parts[i] if i < len(parts) else None
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

    # Load BacDive genome metadata
    click.echo("Loading BacDive genome metadata...")
    genomes = []
    with open(bacdive_jsonl) as f:
        for line in f:
            record = json.loads(line)
            # Extract genome accession as primary key
            genome_acc = record.get("Genome seq. accession number", "").strip()
            if genome_acc.startswith("GCA_"):
                # Store all fields from BacDive
                genome_record = {"genome_id": genome_acc}
                # Add all other fields from the record
                for key, value in record.items():
                    # Skip the accession field we already have as genome_id
                    if key != "Genome seq. accession number":
                        # Clean up field name for SQL (replace spaces and special chars)
                        clean_key = key.replace(" ", "_").replace(".", "").replace("(", "").replace(")", "")
                        genome_record[clean_key] = value
                genomes.append(genome_record)

    db["genomes"].insert_all(genomes, pk="genome_id", replace=True)
    click.echo(f"  Inserted {len(genomes)} genomes")

    # Parse Rfam metadata files
    click.echo("Loading Rfam family metadata...")
    families = parse_rfam_table(Path(family_txt))
    db["family"].insert_all(families, pk="rfam_acc", replace=True)
    click.echo(f"  Inserted {len(families)} families")

    click.echo("Loading Rfam clan metadata...")
    clans = parse_rfam_table(Path(clan_txt))
    db["clan"].insert_all(clans, pk="clan_acc", replace=True)
    click.echo(f"  Inserted {len(clans)} clans")

    click.echo("Loading Rfam clan membership...")
    memberships = parse_rfam_table(Path(clan_membership_txt))
    db["clan_membership"].insert_all(memberships, replace=True)
    click.echo(f"  Inserted {len(memberships)} clan memberships")

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
                db["hits"].insert(hit)
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
    db["hits"].create_index(["genome_id"])
    db["hits"].create_index(["target_name"])
    db["hits"].create_index(["query_name"])
    db["hits"].create_index(["e_value"])
    db["hits"].create_index(["score"])
    db["hits"].create_index(["hit_type"])

    # Create foreign key relationships (informational, not enforced by default)
    click.echo("Adding foreign key references...")
    try:
        # hits.target_name -> family.rfam_acc
        db["hits"].add_foreign_key("target_name", "family", "rfam_acc")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->family: {e}", err=True)

    try:
        # hits.genome_id -> genomes.genome_id
        db["hits"].add_foreign_key("genome_id", "genomes", "genome_id")
    except Exception as e:
        click.echo(f"  Note: Could not add foreign key for hits->genomes: {e}", err=True)

    # Create empty alignment tables (will be populated by create_alignments.py)
    click.echo("Creating alignment table schemas...")

    # Create alignments table
    db["alignments"].create({
        "alignment_id": int,
        "family_id": str,
        "genome_id": str,
        "hit_id": str,  # Reference to hit (genome_id:query_name:seq_from-seq_to)
        "aligned_sequence": str,
        "alignment_start": int,
        "alignment_end": int,
    }, pk="alignment_id", if_not_exists=True)

    # Create alignment_metadata table
    db["alignment_metadata"].create({
        "family_id": str,
        "num_sequences": int,
        "alignment_length": int,
        "consensus_structure": str,
        "created_date": str,
    }, pk="family_id", if_not_exists=True)

    click.echo("  Created alignments and alignment_metadata tables (empty)")

    click.echo(f"\nDatabase created successfully: {output_db}")
    click.echo(f"  Tables: {', '.join(db.table_names())}")

    # Print example query
    click.echo("\nExample queries:")
    click.echo("  -- Find riboswitch hits in genomes that grow below 20°C:")
    click.echo("  SELECT h.*, g.* FROM hits h JOIN genomes g ON h.genome_id = g.genome_id")
    click.echo("  WHERE h.target_name LIKE '%riboswitch%' AND g.Temperature < 20;")
    click.echo("\n  -- Count hits by genome:")
    click.echo("  SELECT g.genome_id, COUNT(*) as hit_count FROM hits h")
    click.echo("  JOIN genomes g ON h.genome_id = g.genome_id GROUP BY g.genome_id;")
    click.echo("\n  -- View alignment for a specific family:")
    click.echo("  SELECT * FROM alignments WHERE family_id = 'RF00001';")


if __name__ == '__main__':
    main()
