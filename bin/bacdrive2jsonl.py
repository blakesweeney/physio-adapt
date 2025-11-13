#!/usr/bin/env python3
import csv
import json
import re
import click
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple


def parse_temperature(temp_str: str) -> Optional[Tuple[float, float]]:
    """
    Parse temperature string into (min, max) tuple.

    Examples:
        '20' -> (20.0, 20.0)
        '25-41' -> (25.0, 41.0)
        '<10' -> (None, 10.0)
        '>45' -> (45.0, None)

    Returns:
        Tuple of (min, max) or None if unparseable
    """
    if not temp_str or not isinstance(temp_str, str):
        return None

    temp_str = temp_str.strip()

    # Handle range format: "25-41"
    range_match = re.match(r'^(\d+(?:\.\d+)?)\s*-\s*(\d+(?:\.\d+)?)$', temp_str)
    if range_match:
        min_val = float(range_match.group(1))
        max_val = float(range_match.group(2))
        return (min_val, max_val)

    # Handle less than: "<10"
    less_than_match = re.match(r'^<\s*(\d+(?:\.\d+)?)$', temp_str)
    if less_than_match:
        return (None, float(less_than_match.group(1)))

    # Handle greater than: ">45"
    greater_than_match = re.match(r'^>\s*(\d+(?:\.\d+)?)$', temp_str)
    if greater_than_match:
        return (float(greater_than_match.group(1)), None)

    # Handle single value: "20"
    single_match = re.match(r'^(\d+(?:\.\d+)?)$', temp_str)
    if single_match:
        val = float(single_match.group(1))
        return (val, val)

    # Could not parse
    return None


def parse_bacdive(file_path: Path) -> List[Dict[str, Any]]:
    """
    Parse a BacDive CSV where empty leading columns inherit values from the previous non-empty row.

    Transforms temperature data into structured format:
    - Parses temperature strings into min/max values
    - Separates genome, strain, and measurement data

    Args:
        file_path: Path to the CSV file

    Returns:
        List of dictionaries containing the parsed and structured data
    """
    results = []
    previous_values = {}

    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames

        for row in reader:
            # Create a new row with inherited values
            raw_row = {}

            for header in headers:
                value = row[header].strip() if row[header] else ""

                if value:
                    # If there's a value, use it and update the cache
                    raw_row[header] = value
                    previous_values[header] = value
                else:
                    # If empty, inherit from previous row if available
                    if header in previous_values:
                        raw_row[header] = previous_values[header]
                    else:
                        raw_row[header] = ""

            # Now restructure the data into genome, strain, and measurement components
            structured_row = {
                "genome": {
                    "accession": raw_row.get("Genome seq. accession number", ""),
                    "database": raw_row.get("Genome Sequence database", "")
                },
                "strain": {
                    "bacdive_id": raw_row.get("ID", ""),
                    "designation": raw_row.get("designation_header", ""),
                    "strain_numbers": raw_row.get("strain_number_header", ""),
                    "is_type_strain": raw_row.get("is_type_strain_header", "") == "1",
                    "species": raw_row.get("species", "")
                },
                "measurement": None
            }

            # Parse temperature measurement if present
            temp_str = raw_row.get("Temperature", "")
            if temp_str:
                parsed_temp = parse_temperature(temp_str)
                if parsed_temp:
                    temp_min, temp_max = parsed_temp
                    structured_row["measurement"] = {
                        "property_type": "temperature",
                        "measurement_kind": raw_row.get("Kind of temperature", ""),
                        "value_min": temp_min,
                        "value_max": temp_max,
                        "value_unit": "°C",
                        "test_result": raw_row.get("Testresult (temperature)", "")
                    }

            results.append(structured_row)

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