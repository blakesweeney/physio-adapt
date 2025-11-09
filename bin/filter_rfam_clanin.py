#!/usr/bin/env python3
"""
Filter Rfam clanin file to only include clans with selected families.
"""

import click
from pathlib import Path


@click.command()
@click.argument("clanin_file", type=click.Path(exists=True, path_type=Path))
@click.argument("families_file", type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
def filter_clanin(clanin_file: Path, families_file: Path, output_file: Path):
    """
    Filter Rfam clanin file to only include clans with selected families.

    The clanin file format has lines like:
    CL00001 RF00001 RF00002 RF00003 ...

    Where the first column is the clan ID and subsequent columns are family accessions.
    We keep only clans that contain at least one of our selected families.

    CLANIN_FILE: Path to Rfam.clanin file
    FAMILIES_FILE: File with selected family accessions (one per line)
    OUTPUT_FILE: Path to output filtered clanin file
    """
    # Read selected families into a set for fast lookup
    selected_families = set()
    with open(families_file, "r") as f:
        for line in f:
            family = line.strip()
            if family:
                selected_families.add(family)

    click.echo(f"Loaded {len(selected_families)} selected families")

    # Filter clanin file
    clans_kept = 0
    clans_total = 0

    with open(clanin_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            clans_total += 1
            parts = line.strip().split()

            if len(parts) < 2:
                continue

            _clan_id = parts[0]
            families_in_clan = set(parts[1:])

            # Check if this clan contains any of our selected families
            if families_in_clan & selected_families:
                outfile.write(line)
                clans_kept += 1

    click.echo(f"Kept {clans_kept} clans out of {clans_total} total clans")
    click.echo(f"Filtered clanin file written to {output_file}")


if __name__ == "__main__":
    filter_clanin()
