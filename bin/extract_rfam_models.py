#!/usr/bin/env python3
"""
Extract Rfam covariance models by RNA type.
"""

import gzip
import subprocess
import sys
import click
from pathlib import Path


def parse_family_file(
    family_gz_path: Path, rna_types: set[str]
) -> dict[str, list[str]]:
    """
    Parse Rfam family.txt.gz and extract accessions for specified RNA types.

    Returns a dictionary mapping RNA type to list of family accessions.
    """
    type_to_families = {rna_type: [] for rna_type in rna_types}
    all_types = set()

    with gzip.open(family_gz_path, "rt") as f:
        # Skip header line
        _header = next(f)

        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 18:
                continue

            rfam_acc = parts[0]  # Column 1: rfam_acc
            rna_type = parts[17]  # Column 18: type (0-indexed = 17)

            all_types.add(rna_type)

            if rna_type in rna_types:
                type_to_families[rna_type].append(rfam_acc)

    # Check if we found any families
    for rna_type in rna_types:
        if not type_to_families[rna_type]:
            click.echo(
                f"WARNING: No families found for RNA type '{rna_type}'", err=True
            )
            click.echo(f"Available types: {', '.join(sorted(all_types))}", err=True)

    return type_to_families


@click.command()
@click.argument("rfam_cm_gz", type=click.Path(exists=True, path_type=Path))
@click.argument("family_txt_gz", type=click.Path(exists=True, path_type=Path))
@click.argument("rna_types_file", type=click.Path(exists=True, path_type=Path))
@click.argument("output_cm", type=click.Path(path_type=Path))
@click.argument("output_families", type=click.Path(path_type=Path))
def extract_models(
    rfam_cm_gz: Path,
    family_txt_gz: Path,
    rna_types_file: Path,
    output_cm: Path,
    output_families: Path,
):
    """
    Extract Rfam covariance models for specified RNA types.

    RFAM_CM_GZ: Path to Rfam.cm.gz file
    FAMILY_TXT_GZ: Path to Rfam family.txt.gz file
    RNA_TYPES_FILE: File with RNA types (one per line)
    OUTPUT_CM: Path to output CM file (will be compressed with cmpress)
    OUTPUT_FAMILIES: Path to output file listing extracted families
    """
    # Read RNA types from file
    with open(rna_types_file, "r") as f:
        rna_types = {
            line.strip() for line in f if line.strip() and not line.startswith("#")
        }

    if not rna_types:
        click.echo("ERROR: No RNA types specified in input file", err=True)
        sys.exit(1)

    click.echo(f"Searching for RNA types: {', '.join(sorted(rna_types))}")

    # Parse family file to get accessions
    type_to_families = parse_family_file(family_txt_gz, rna_types)

    # Collect all family accessions
    all_families = []
    for rna_type, families in type_to_families.items():
        click.echo(f"Found {len(families)} families for type '{rna_type}'")
        all_families.extend(families)

    if not all_families:
        click.echo("ERROR: No families found for any specified RNA type", err=True)
        sys.exit(1)

    # Write families list to output file
    with open(output_families, "w") as f:
        for family in sorted(all_families):
            f.write(family + "\n")

    click.echo(f"Extracting {len(all_families)} models from Rfam.cm.gz...")

    # Decompress Rfam.cm.gz
    click.echo("Decompressing Rfam.cm.gz...")
    subprocess.run(
        ["gunzip", "-c", str(rfam_cm_gz)], stdout=open("Rfam.cm", "w"), check=True
    )

    # Extract models using cmfetch
    click.echo("Extracting models with cmfetch...")
    with open(output_cm, "w") as out:
        for family_acc in all_families:
            try:
                result = subprocess.run(
                    ["cmfetch", "Rfam.cm", family_acc],
                    capture_output=True,
                    text=True,
                    check=True,
                )
                out.write(result.stdout)
            except subprocess.CalledProcessError as e:
                click.echo(f"WARNING: Could not fetch {family_acc}: {e}", err=True)

    # Press the CM file for faster searching
    click.echo("Pressing CM file with cmpress...")
    subprocess.run(["cmpress", "-F", str(output_cm)], check=True)

    click.echo(f"Successfully created {output_cm} with {len(all_families)} models")


if __name__ == "__main__":
    extract_models()
