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
    query = "SELECT DISTINCT target_name FROM hits ORDER BY target_name"
    return [row['target_name'] for row in db.execute(query).fetchall()]


def get_hits_for_family(db: Database, family_id: str) -> List[Dict]:
    """Get all hits for a specific family."""
    query = """
    SELECT
        h.genome_id,
        h.query_name,
        h.seq_from,
        h.seq_to,
        h.hit_sequence,
        h.target_accession,
        h.strand
    FROM hits h
    WHERE h.target_name = ?
    AND h.hit_sequence IS NOT NULL
    ORDER BY h.genome_id, h.query_name, h.seq_from
    """
    return list(db.execute(query, [family_id]).fetchall())


def create_fasta_from_hits(hits: List[Dict], output_path: Path) -> None:
    """Create FASTA file from hit sequences."""
    with open(output_path, 'w') as f:
        for hit in hits:
            # Create unique ID: genome_id:query_name:seq_from-seq_to
            hit_id = f"{hit['genome_id']}:{hit['query_name']}:{hit['seq_from']}-{hit['seq_to']}"
            f.write(f">{hit_id}\n")
            f.write(f"{hit['hit_sequence']}\n")


def run_cmalign(cm_file: Path, fasta_file: Path, output_sto: Path, family_id: str) -> bool:
    """Run cmalign to create alignment."""
    try:
        # First, extract the CM model for this family
        cm_extract = output_sto.parent / f"{family_id}.cm"

        click.echo(f"    Extracting CM model for {family_id}...")
        subprocess.run(
            ['cmfetch', str(cm_file), family_id],
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
    family_id: str,
    aligned_seqs: Dict[str, str],
    consensus_structure: str,
    alignment_length: int
) -> None:
    """Store alignment data in database."""

    # Store alignment metadata
    metadata = {
        'family_id': family_id,
        'num_sequences': len(aligned_seqs),
        'alignment_length': alignment_length,
        'consensus_structure': consensus_structure,
        'created_date': datetime.now().isoformat()
    }
    db['alignment_metadata'].insert(metadata, replace=True)

    # Store individual aligned sequences
    alignments = []
    for seq_id, aligned_seq in aligned_seqs.items():
        # Parse hit_id: genome_id:query_name:seq_from-seq_to
        parts = seq_id.split(':')
        if len(parts) >= 3:
            genome_id = parts[0]

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
                'family_id': family_id,
                'genome_id': genome_id,
                'hit_id': seq_id,
                'aligned_sequence': aligned_seq,
                'alignment_start': align_start,
                'alignment_end': align_end
            })

    # Batch insert all alignments for this family
    db['alignments'].insert_all(alignments, alter=True)
    click.echo(f"    Stored {len(alignments)} aligned sequences")


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

        for i, family_id in enumerate(families, 1):
            click.echo(f"\n[{i}/{len(families)}] Processing {family_id}...")

            # Get hits for this family
            hits = get_hits_for_family(db, family_id)
            click.echo(f"  Found {len(hits)} hits")

            if len(hits) == 0:
                click.echo(f"  Skipping {family_id} (no sequences)")
                continue

            # Create FASTA file in tmpdir
            fasta_file = tmpdir / f"{family_id}.fa"
            create_fasta_from_hits(hits, fasta_file)

            # Optionally save FASTA file to output directory
            if save_fasta:
                output_fasta = output_path / f"{family_id}.fa"
                import shutil
                shutil.copy2(fasta_file, output_fasta)
                click.echo(f"  Saved FASTA: {output_fasta.name}")

            # Run cmalign - save Stockholm to output directory
            sto_file = output_path / f"{family_id}.sto"
            if run_cmalign(Path(rfam_cm), fasta_file, sto_file, family_id):
                # Parse alignment
                try:
                    aligned_seqs, consensus_structure, alignment_length = parse_stockholm_alignment(sto_file)

                    # Store in database
                    store_alignment_in_db(db, family_id, aligned_seqs, consensus_structure, alignment_length)

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
                click.echo(f"  ✗ cmalign failed for {family_id}")

    click.echo(f"\n{'='*60}")
    click.echo("Alignment generation complete!")
    click.echo(f"  Successful: {success_count}/{len(families)}")
    click.echo(f"  Errors: {error_count}/{len(families)}")
    click.echo(f"\nDatabase updated: {database}")


if __name__ == '__main__':
    main()
