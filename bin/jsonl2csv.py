#!/usr/bin/env python3
"""
Convert BacDive JSONL format to flattened CSV for R/tidyverse analysis.

Usage:
    python bin/jsonl2csv.py input.jsonl output.csv
"""

import csv
import json
from pathlib import Path

import click


def flatten_record(record):
    """
    Flatten nested JSON record into a single-level dictionary.

    Args:
        record: Dict with nested structure (genome, strain, measurement)

    Returns:
        Dict with flattened keys like 'genome_accession', 'strain_species', etc.
        None/null values are converted to empty strings for CSV compatibility.
    """
    flattened = {}

    # Flatten genome fields
    if "genome" in record and record["genome"] is not None:
        for key, value in record["genome"].items():
            flattened[f"genome_{key}"] = "" if value is None else value

    # Flatten strain fields
    if "strain" in record and record["strain"] is not None:
        for key, value in record["strain"].items():
            flattened[f"strain_{key}"] = "" if value is None else value

    # Flatten measurement fields
    if "measurement" in record and record["measurement"] is not None:
        for key, value in record["measurement"].items():
            flattened[f"measurement_{key}"] = "" if value is None else value

    return flattened


def get_all_fieldnames(jsonl_path):
    """
    Scan the JSONL file to discover all possible field names.
    This ensures we have all columns even if some records have missing fields.

    Args:
        jsonl_path: Path to input JSONL file

    Returns:
        Sorted list of all unique field names
    """
    all_fields = set()

    with open(jsonl_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                record = json.loads(line)
                flattened = flatten_record(record)
                all_fields.update(flattened.keys())

    # Sort fields for consistent column ordering
    # Group by prefix (genome_, strain_, measurement_) for readability
    genome_fields = sorted([f for f in all_fields if f.startswith("genome_")])
    strain_fields = sorted([f for f in all_fields if f.startswith("strain_")])
    measurement_fields = sorted([f for f in all_fields if f.startswith("measurement_")])

    return genome_fields + strain_fields + measurement_fields


@click.command()
@click.argument("input_jsonl", type=click.Path(exists=True, path_type=Path))
@click.argument("output_csv", type=click.Path(path_type=Path))
def jsonl2csv(input_jsonl, output_csv):
    """
    Convert BacDive JSONL to flattened CSV format.

    Reads INPUT_JSONL and writes flattened CSV to OUTPUT_CSV.
    All nested fields are flattened with underscore notation
    (e.g., genome.accession becomes genome_accession).
    """
    # First pass: discover all fields
    click.echo(f"Scanning {input_jsonl} to discover fields...")
    fieldnames = get_all_fieldnames(input_jsonl)
    click.echo(f"Found {len(fieldnames)} unique fields")

    # Second pass: write CSV
    click.echo(f"Converting to CSV: {output_csv}")
    records_written = 0

    with open(input_jsonl, "r", encoding="utf-8") as infile, \
         open(output_csv, "w", encoding="utf-8", newline="") as outfile:

        writer = csv.DictWriter(outfile, fieldnames=fieldnames,
                               extrasaction="ignore", restval="")
        writer.writeheader()

        for line in infile:
            if line.strip():
                record = json.loads(line)
                flattened = flatten_record(record)
                writer.writerow(flattened)
                records_written += 1

    click.echo(f"Successfully wrote {records_written} records to {output_csv}")


if __name__ == "__main__":
    jsonl2csv()
