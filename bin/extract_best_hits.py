#!/usr/bin/env python3
"""
Extract best scoring hit per Rfam model from cmscan tblout file.
"""
import click
from pathlib import Path


def parse_tblout_line(line: str) -> dict:
    """
    Parse a single line from cmscan tblout output (--fmt 2).

    Format columns:
    1. target name
    2. target accession
    3. query name
    4. query accession
    5. mdl (model type)
    6. mdl from
    7. mdl to
    8. seq from
    9. seq to
    10. strand
    11. trunc
    12. pass
    13. gc
    14. bias
    15. score
    16. E-value
    17. inc
    18. description (rest of line)
    """
    parts = line.strip().split(None, 17)  # Split on whitespace, max 18 parts

    if len(parts) < 17:
        return None

    return {
        'target_name': parts[0],
        'target_acc': parts[1],
        'query_name': parts[2],
        'query_acc': parts[3],
        'mdl': parts[4],
        'mdl_from': parts[5],
        'mdl_to': parts[6],
        'seq_from': parts[7],
        'seq_to': parts[8],
        'strand': parts[9],
        'trunc': parts[10],
        'pass': parts[11],
        'gc': parts[12],
        'bias': parts[13],
        'score': float(parts[14]),
        'e_value': parts[15],
        'inc': parts[16],
        'description': parts[17] if len(parts) > 17 else '',
        'original_line': line
    }


@click.command()
@click.argument('tblout_file', type=click.Path(exists=True, path_type=Path))
@click.argument('output_file', type=click.Path(path_type=Path))
@click.option('--genome-id', default='genome', help='Genome identifier for output')
def extract_best_hits(tblout_file: Path, output_file: Path, genome_id: str):
    """
    Extract best scoring hit per Rfam model from cmscan tblout file.

    Finds the highest scoring hit for each model (target), even if it's below
    the reporting threshold.

    TBLOUT_FILE: Path to cmscan tblout output file
    OUTPUT_FILE: Path to output file with best hits
    """
    # Dictionary to store best hit for each model
    # Key: target_acc (model accession), Value: parsed hit dictionary
    best_hits = {}

    total_hits = 0
    with open(tblout_file, 'r') as f:
        for line in f:
            # Skip comment and header lines
            if line.startswith('#'):
                continue

            hit = parse_tblout_line(line)
            if hit is None:
                continue

            total_hits += 1
            target_acc = hit['target_acc']
            score = hit['score']

            # Keep the best scoring hit for each model
            if target_acc not in best_hits or score > best_hits[target_acc]['score']:
                best_hits[target_acc] = hit

    # Write output
    with open(output_file, 'w') as out:
        # Write header
        out.write(f"# Best scoring hit per Rfam model for genome: {genome_id}\n")
        out.write(f"# Total hits processed: {total_hits}\n")
        out.write(f"# Unique models with hits: {len(best_hits)}\n")
        out.write("#\n")
        out.write("# Columns: target_name, target_acc, query_name, query_acc, mdl, mdl_from, mdl_to, ")
        out.write("seq_from, seq_to, strand, trunc, pass, gc, bias, score, e_value, inc, description\n")
        out.write("#\n")

        # Sort by target accession and write
        for target_acc in sorted(best_hits.keys()):
            out.write(best_hits[target_acc]['original_line'])

    click.echo(f"Processed {total_hits} total hits")
    click.echo(f"Found best hits for {len(best_hits)} unique models")
    click.echo(f"Output written to {output_file}")


if __name__ == '__main__':
    extract_best_hits()
