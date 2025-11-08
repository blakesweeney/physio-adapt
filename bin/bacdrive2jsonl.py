#!/usr/bin/env python3
import csv
import json
import click
from pathlib import Path
from typing import List, Dict, Any


def parse_bacdive(file_path: Path) -> List[Dict[str, Any]]:
    """
    Parse a BacDive CSV where empty leading columns inherit values from the previous non-empty row.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        List of dictionaries containing the parsed data
    """
    results = []
    previous_values = {}
    
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        
        for row in reader:
            # Create a new row with inherited values
            processed_row = {}
            
            for header in headers:
                value = row[header].strip() if row[header] else ""
                
                if value:
                    # If there's a value, use it and update the cache
                    processed_row[header] = value
                    previous_values[header] = value
                else:
                    # If empty, inherit from previous row if available
                    if header in previous_values:
                        processed_row[header] = previous_values[header]
                    else:
                        processed_row[header] = ""
            
            results.append(processed_row)
    
    return results


@click.command()
@click.argument('filename', type=click.Path(exists=True, path_type=Path))
@click.argument('output', type=click.Path(path_type=Path))
def bacdive2json(filename: Path, output: Path):
    """
    Convert BacDive CSV to JSONL format.
    
    FILENAME: Path to the input CSV file
    OUTPUT: Path to the output JSONL file
    """
    try:
        # Parse the CSV file
        records = parse_bacdive(filename)
        
        # Write to JSONL
        with open(output, 'w', encoding='utf-8') as f:
            for record in records:
                json.dump(record, f, ensure_ascii=False)
                f.write('\n')
        
        click.echo(f"Successfully converted {len(records)} records from {filename} to {output}")
        
    except Exception as e:
        click.echo(f"Error: {e}", err=True)
        raise click.Abort()


if __name__ == '__main__':
    bacdive2json()