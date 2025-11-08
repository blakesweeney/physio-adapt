#!/usr/bin/env python3
"""
Extract GCA accession numbers from BacDive JSONL file.
"""
import json
import click
from pathlib import Path


@click.command()
@click.argument('jsonl_file', type=click.Path(exists=True, path_type=Path))
@click.argument('output', type=click.Path(path_type=Path))
@click.option('--limit', type=int, default=None, help='Limit number of accessions')
def extract_accessions(jsonl_file: Path, output: Path, limit: int):
    """
    Extract GCA accession numbers from JSONL file.

    JSONL_FILE: Path to the input JSONL file
    OUTPUT: Path to output text file with one accession per line
    """
    accessions = set()

    with open(jsonl_file, 'r') as f:
        for line in f:
            record = json.loads(line)
            # Extract GCA accessions from "Genome seq. accession number" field
            acc = record.get("Genome seq. accession number", "").strip()
            if acc.startswith("GCA_"):
                accessions.add(acc)

    # Sort for reproducibility
    accessions_list = sorted(accessions)

    # Apply limit if specified
    if limit:
        accessions_list = accessions_list[:limit]

    with open(output, 'w') as f:
        for acc in accessions_list:
            f.write(acc + "\n")

    click.echo(f"Found {len(accessions)} unique GCA accessions, wrote {len(accessions_list)} to {output}")


if __name__ == '__main__':
    extract_accessions()
