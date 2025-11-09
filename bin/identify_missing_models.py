#!/usr/bin/env python3
"""
Identify which models from a selected set have no hits in a tblout file.

This script compares the selected families (RNA type subset) against
the models that have hits in a cmscan tblout file, and outputs the
list of models that are missing (have zero hits).
"""

import click
from pathlib import Path


def parse_tblout_models(tblout_path: Path) -> set[str]:
    """
    Parse tblout file and return set of unique model accessions that have hits.

    Returns set of target accessions (column 2 in tblout format).
    """
    models_with_hits = set()

    with open(tblout_path) as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue

            parts = line.strip().split()
            if len(parts) < 2:
                continue

            # Column 2 is target accession (0-indexed = 1)
            target_acc = parts[1]
            models_with_hits.add(target_acc)

    return models_with_hits


def read_selected_families(families_path: Path) -> set[str]:
    """
    Read selected families file (one accession per line).

    Returns set of family accessions.
    """
    families = set()

    with open(families_path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                families.add(line)

    return families


@click.command()
@click.argument('tblout_file', type=click.Path(exists=True, path_type=Path))
@click.argument('selected_families_file', type=click.Path(exists=True, path_type=Path))
@click.argument('output_file', type=click.Path(path_type=Path))
@click.option('--genome-id', default='genome', help='Genome identifier for logging')
def main(tblout_file: Path, selected_families_file: Path, output_file: Path, genome_id: str):
    """
    Identify models from selected families that have no hits in tblout.

    TBLOUT_FILE: cmscan tblout output file
    SELECTED_FAMILIES_FILE: File with selected family accessions (one per line)
    OUTPUT_FILE: Output file with missing model accessions (one per line)
    """
    # Read selected families
    selected_families = read_selected_families(selected_families_file)
    click.echo(f"Selected families: {len(selected_families)}")

    # Parse tblout to find models with hits
    models_with_hits = parse_tblout_models(tblout_file)
    click.echo(f"Models with hits in {genome_id}: {len(models_with_hits)}")

    # Find missing models (selected but not in hits)
    missing_models = selected_families - models_with_hits
    click.echo(f"Missing models: {len(missing_models)}")

    # Write missing models to output
    with open(output_file, 'w') as out:
        for model_acc in sorted(missing_models):
            out.write(f"{model_acc}\n")

    if missing_models:
        click.echo(f"Missing models written to {output_file}")
    else:
        click.echo(f"All selected models have hits - created empty file: {output_file}")


if __name__ == '__main__':
    main()
