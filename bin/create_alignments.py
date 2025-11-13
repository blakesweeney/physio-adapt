#!/usr/bin/env python3
"""
Create multiple sequence alignments for all Rfam families with hits.

This script:
1. Queries the database to find all families with hits
2. For each family, extracts hit sequences
3. Runs cmalign to create Stockholm format alignment
4. Parses the alignment and stores in database
5. Stores alignment metadata
"""

import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import click
from Bio import AlignIO
from sqlite_utils import Database


def get_families_with_hits(db: Database) -> List[str]:
    """Query database to find all Rfam families that have hits."""
    query = "SELECT DISTINCT rfam_acc FROM hits ORDER BY rfam_acc"
    return [row[0] for row in db.execute(query).fetchall()]


def get_hits_for_family(db: Database, rfam_acc: str) -> List[Dict]:
    """Get all hits for a specific family."""
    query = """
    SELECT
        h.genome_id,
        h.query_name,
        h.seq_from,
        h.seq_to,
        h.hit_sequence,
        h.strand
    FROM hits h
    WHERE h.rfam_acc = ?
    AND h.hit_sequence IS NOT NULL
    ORDER BY h.genome_id, h.query_name, h.seq_from
    """
    rows = db.execute(query, [rfam_acc]).fetchall()
    # Convert tuples to dictionaries
    keys = ['genome_id', 'query_name', 'seq_from', 'seq_to', 'hit_sequence', 'strand']
    return [dict(zip(keys, row)) for row in rows]


def create_fasta_from_hits(hits: List[Dict], output_path: Path) -> None:
    """Create FASTA file from hit sequences.

    Sequence IDs are in format: accession/start-stop
    """
    with open(output_path, 'w') as f:
        for hit in hits:
            # Create unique ID: query_name/seq_from-seq_to
            hit_id = f"{hit['query_name']}/{hit['seq_from']}-{hit['seq_to']}"
            f.write(f">{hit_id}\n")
            f.write(f"{hit['hit_sequence']}\n")


def run_cmalign(cm_file: Path, fasta_file: Path, output_sto: Path, rfam_acc: str) -> bool:
    """Run cmalign to create alignment."""
    try:
        # First, extract the CM model for this family
        cm_extract = output_sto.parent / f"{rfam_acc}.cm"

        click.echo(f"    Extracting CM model for {rfam_acc}...")
        subprocess.run(
            ['cmfetch', str(cm_file), rfam_acc],
            stdout=open(cm_extract, 'w'),
            stderr=subprocess.PIPE,
            check=True
        )

        # Run cmalign
        click.echo("    Running cmalign...")
        subprocess.run(
            ['cmalign', '-o', str(output_sto), str(cm_extract), str(fasta_file)],
            stderr=subprocess.PIPE,
            check=True
        )

        # Clean up extracted CM
        cm_extract.unlink()

        return True

    except subprocess.CalledProcessError as e:
        click.echo(f"    Error running cmalign: {e.stderr.decode()}", err=True)
        return False


def parse_stockholm_alignment(sto_file: Path) -> Tuple[Dict, str, int]:
    """
    Parse Stockholm alignment file.

    Returns:
        - Dictionary mapping sequence IDs to aligned sequences
        - Consensus structure (if available)
        - Alignment length
    """
    alignment = AlignIO.read(sto_file, 'stockholm')

    aligned_seqs = {}
    for record in alignment:
        aligned_seqs[record.id] = str(record.seq)

    # Try to get consensus structure from Stockholm annotation
    consensus_structure = ""
    if hasattr(alignment, 'column_annotations') and 'secondary_structure' in alignment.column_annotations:
        consensus_structure = alignment.column_annotations['secondary_structure']

    alignment_length = alignment.get_alignment_length()

    return aligned_seqs, consensus_structure, alignment_length


def store_alignment_in_db(
    db: Database,
    rfam_acc: str,
    aligned_seqs: Dict[str, str],
    consensus_structure: str,
    alignment_length: int,
    sto_file: Path,
    genome_lookup: Dict[str, str]
) -> None:
    """Store alignment data in database.

    Args:
        db: Database connection
        rfam_acc: Rfam family accession
        aligned_seqs: Dictionary of sequence IDs to aligned sequences
        consensus_structure: Consensus secondary structure
        alignment_length: Length of alignment
        sto_file: Path to Stockholm file (to store as blob)
        genome_lookup: Dictionary mapping accessions to genome IDs
    """

    # Store alignment metadata
    metadata = {
        'rfam_acc': rfam_acc,
        'num_sequences': len(aligned_seqs),
        'alignment_length': alignment_length,
        'consensus_structure': consensus_structure,
        'created_date': datetime.now().isoformat()
    }
    db['alignment_metadata'].insert(metadata, replace=True, strict=True)

    # Store complete Stockholm file as blob
    with open(sto_file, 'rb') as f:
        stockholm_blob = f.read()

    alignment_file_record = {
        'rfam_acc': rfam_acc,
        'stockholm_alignment': stockholm_blob,
        'file_size': len(stockholm_blob),
        'created_date': datetime.now().isoformat()
    }
    db['alignment_files'].insert(alignment_file_record, replace=True, strict=True)

    # Store individual aligned sequences
    alignments = []
    for seq_id, aligned_seq in aligned_seqs.items():
        # Parse hit_id: accession/seq_from-seq_to
        # Extract accession to look up genome_id
        accession = seq_id.split('/')[0] if '/' in seq_id else seq_id
        genome_id = genome_lookup.get(accession, None)

        # Get alignment coordinates (1-indexed, without gaps)
        # Find first and last non-gap positions
        ungapped_pos = [i for i, c in enumerate(aligned_seq) if c != '-' and c != '.']
        if ungapped_pos:
            align_start = ungapped_pos[0] + 1
            align_end = ungapped_pos[-1] + 1
        else:
            align_start = 0
            align_end = 0

        alignments.append({
            'rfam_acc': rfam_acc,
            'genome_id': genome_id,
            'hit_id': seq_id,
            'aligned_sequence': aligned_seq,
            'alignment_start': align_start,
            'alignment_end': align_end
        })

    # Batch insert all alignments for this family
    db['alignments'].insert_all(alignments, alter=True, strict=True)
    click.echo(f"    Stored {len(alignments)} aligned sequences and Stockholm file ({len(stockholm_blob)} bytes)")


@click.command()
@click.argument('database', type=click.Path(exists=True))
@click.argument('rfam_cm', type=click.Path(exists=True))
@click.option('--output-dir', type=click.Path(), default='.', help='Directory to save Stockholm alignment files')
@click.option('--save-fasta', is_flag=True, help='Also save FASTA files used for alignment')
@click.option('--limit', type=int, default=None, help='Limit number of families to process (for testing)')
def main(database: str, rfam_cm: str, output_dir: str, save_fasta: bool, limit: int):
    """
    Create alignments for all Rfam families with hits in the database.

    DATABASE: Path to SQLite database file
    RFAM_CM: Path to Rfam.cm file (must have cmfetch index)
    """
    db = Database(database)

    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    click.echo(f"Saving Stockholm files to: {output_path.absolute()}")

    # Build genome lookup from sequences table
    click.echo("Building accession to genome lookup...")
    genome_lookup = {}
    for row in db.execute("SELECT accession, genome_id FROM sequences").fetchall():
        genome_lookup[row[0]] = row[1]
    click.echo(f"  Loaded {len(genome_lookup)} sequence-to-genome mappings")

    # Get all families with hits
    click.echo("Finding families with hits...")
    families = get_families_with_hits(db)
    click.echo(f"Found {len(families)} families with hits")

    if limit:
        families = families[:limit]
        click.echo(f"Limiting to {limit} families for testing")

    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Process each family
        success_count = 0
        error_count = 0

        for i, rfam_acc in enumerate(families, 1):
            click.echo(f"\n[{i}/{len(families)}] Processing {rfam_acc}...")

            # Get hits for this family
            hits = get_hits_for_family(db, rfam_acc)
            click.echo(f"  Found {len(hits)} hits")

            if len(hits) == 0:
                click.echo(f"  Skipping {rfam_acc} (no sequences)")
                continue

            # Create FASTA file in tmpdir
            fasta_file = tmpdir / f"{rfam_acc}.fa"
            create_fasta_from_hits(hits, fasta_file)

            # Optionally save FASTA file to output directory
            if save_fasta:
                output_fasta = output_path / f"{rfam_acc}.fa"
                import shutil
                shutil.copy2(fasta_file, output_fasta)
                click.echo(f"  Saved FASTA: {output_fasta.name}")

            # Run cmalign - save Stockholm to output directory
            sto_file = output_path / f"{rfam_acc}.sto"
            if run_cmalign(Path(rfam_cm), fasta_file, sto_file, rfam_acc):
                # Parse alignment
                try:
                    aligned_seqs, consensus_structure, alignment_length = parse_stockholm_alignment(sto_file)

                    # Store in database
                    store_alignment_in_db(db, rfam_acc, aligned_seqs, consensus_structure, alignment_length, sto_file, genome_lookup)

                    success_count += 1
                    click.echo(f"  ✓ Successfully aligned {len(aligned_seqs)} sequences (length: {alignment_length})")
                    click.echo(f"  Saved Stockholm: {sto_file.name}")

                except Exception as e:
                    click.echo(f"  ✗ Error parsing/storing alignment: {e}", err=True)
                    error_count += 1
                    # Remove failed Stockholm file
                    if sto_file.exists():
                        sto_file.unlink()
            else:
                error_count += 1
                click.echo(f"  ✗ cmalign failed for {rfam_acc}")

    click.echo(f"\n{'='*60}")
    click.echo("Alignment generation complete!")
    click.echo(f"  Successful: {success_count}/{len(families)}")
    click.echo(f"  Errors: {error_count}/{len(families)}")

    # Add foreign key constraints for alignment tables
    click.echo("\nAdding foreign key constraints for alignment tables...")
    try:
        db["alignments"].add_foreign_key("rfam_acc", "family", "rfam_acc", ignore=True)
        db["alignments"].add_foreign_key("genome_id", "genomes", "genome_id", ignore=True)
        db["alignment_metadata"].add_foreign_key("rfam_acc", "family", "rfam_acc", ignore=True)
        db["alignment_files"].add_foreign_key("rfam_acc", "family", "rfam_acc", ignore=True)
        click.echo("  Foreign keys added successfully")
    except Exception as e:
        click.echo(f"  Warning: Could not add all foreign keys: {e}", err=True)

    click.echo(f"\nDatabase updated: {database}")


if __name__ == '__main__':
    main()
