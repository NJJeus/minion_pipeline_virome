import argparse
import os
import re
from pathlib import Path
from collections import defaultdict

# Note: This script requires pandas and Biopython to be installed:
# pip install pandas biopython

# Using try/except for mandatory external library imports
try:
    import pandas as pd
    from Bio import SeqIO
except ImportError:
    print("Error: This script requires the 'pandas' and 'biopython' libraries.")
    print("Please install them using: pip install pandas biopython")
    exit(1)

def get_sample_name(filename):
    """
    Extracts the sample name from a filename based on the expected pattern.
    Assumes the sample name is everything before the first underscore or the file extension.
    """
    stem = Path(filename).stem
    # Handle both {sample}.csv and {sample}_matched_contigs.fasta
    match = re.match(r'(.+?)(?=_matched_contigs)?$', stem)
    return match.group(1) if match else stem

def clean_taxon_name(taxon_name):
    """Cleans up the Taxon name for use as a filename."""
    # Convert to string, replace spaces with underscores, and strip non-alphanumeric/underscore chars
    clean_name = str(taxon_name).replace(' ', '_').replace('.', '_').strip()
    # Remove any characters that might cause issues in file paths
    clean_name = re.sub(r'[^\w\-_\.]', '', clean_name)
    return clean_name

def aggregate_csv_data(csv_filepath):
    """
    Reads a CSV file, groups by 'Name', and keeps only the row with the 
    maximum 'Score' for each unique 'Name'.
    Returns a dictionary mapping Name to its best metadata.
    """
    print(f"  -> Aggregating data from {Path(csv_filepath).name}...")
    
    try:
        # Read the CSV, ensuring required columns exist
        df = pd.read_csv(csv_filepath)
        required_cols = ['Query', 'Name', 'Taxon', 'Positive terms', 'Score', 'Length_contig']
        
        for col in required_cols:
            if col not in df.columns:
                print(f"  [Error] Missing required column '{col}' in CSV. Skipping file.")
                return None
        
        # Step 2: Group by 'Name' and find the index of the maximum 'Score' within each group.
        # Use .idxmax() to get the index of the row with the max score for each Name.
        max_score_indices = df.groupby('Name')['Score'].idxmax()
        
        # Select the full rows corresponding to those max scores.
        df_max_score = df.loc[max_score_indices, required_cols]
        
        # Create the mapping dictionary
        name_map = {}
        for _, row in df_max_score.iterrows():
            # Store the Taxon name in a clean format for later filename use
            clean_taxon = clean_taxon_name(row['Taxon'])
            
            name_map[row['Name']] = {
                'Taxon': row['Taxon'], # Original Taxon name
                'CleanTaxon': clean_taxon, # Cleaned Taxon name for file grouping
                'Positive terms': row['Positive terms'],
                'Score': row['Score'],
                'Length': row['Length_contig']
            }
            
        print(f"  -> Found {len(name_map)} unique 'Name' entries after aggregation.")
        return name_map

    except Exception as e:
        print(f"  [Error] Failed to read or process CSV file: {e}")
        return None


def process_sample_pair(sample_name, csv_path, fasta_path):
    """
    Processes a single CSV/FASTA pair: aggregates data and modifies FASTA records.
    Returns a list of modified records.
    """
    
    name_map = aggregate_csv_data(csv_path)
    if not name_map:
        return []

    print(f"  -> Processing FASTA file {Path(fasta_path).name}...")
    
    modified_records = []
    
    try:
        # Step 4: Iter through records, matching the ids of the records with Name
        for rec in SeqIO.parse(fasta_path, "fasta"):
            # Use rec.id for matching against the 'Name' column
            rec_name = rec.id.split()[0] # Use only the part before space/description
            
            if rec_name in name_map:
                metadata = name_map[rec_name]
                
                # Use the cleaned name for the new ID construction
                clean_taxon = metadata['CleanTaxon']
                terms = str(metadata['Positive terms']).replace(' ', '_').replace(',', ';')
                score = round(metadata['Score'], 2)
                length = metadata['Length']


                # Step 5: Modify record ID and description
                
                # New ID format: {Sample}_{rec.id}_{CleanTaxon}
                rec.id = f"{sample_name}_{rec_name}_{clean_taxon}"
                
                # New description/header line format
                rec.description = f"Taxon={metadata['Taxon']} Terms={terms} Score={score} Length = {length}"
                
                # Set name to match the new ID
                rec.name = rec.id
                
                # Store the clean Taxon name as an attribute for later grouping
                rec.clean_taxon = clean_taxon
                
                modified_records.append(rec)
            # else: skip records not found in the CSV map
            
        print(f"  -> Modified {len(modified_records)} records for sample '{sample_name}'.")
        return modified_records

    except Exception as e:
        print(f"  [Error] Failed to process FASTA file {Path(fasta_path).name}: {e}")
        return []


def main():
    parser = argparse.ArgumentParser(
        description="Process CSV and FASTA pairs: aggregate data, annotate FASTA records, sort, and output grouped by Taxon."
    )
    parser.add_argument(
        'input_dir', 
        type=Path, 
        help="Path to the input directory containing {sample}.csv and {sample}_matched_contigs.fasta files."
    )
    parser.add_argument(
        'output_dir', 
        type=Path, 
        help="Path to the output directory where annotated and combined FASTA files will be saved."
    )
    
    args = parser.parse_args()
    
    # Setup directories
    input_dir = args.input_dir
    output_dir = args.output_dir
    
    if not input_dir.is_dir():
        print(f"Error: Input directory not found at {input_dir}")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    taxon_output_dir = output_dir / "BY_TAXON"
    taxon_output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing files from: {input_dir}")
    print(f"Outputting results to: {output_dir}")
    print(f"Taxon-grouped output directory: {taxon_output_dir}")
    
    # Use defaultdict to easily group records by Taxon
    taxon_records_map = defaultdict(list)
    
    # 1. Read each CSV table
    for csv_path in input_dir.glob('*.csv'):
        sample_name = get_sample_name(csv_path.name)
        fasta_filename = f"{sample_name}_matched_contigs.fasta"
        fasta_path = input_dir / fasta_filename
        
        print("-" * 40)
        print(f"Found sample: {sample_name}")

        # 3. Open corresponding fasta file
        if fasta_path.is_file():
            
            # Process the pair and get the list of modified records
            modified_records = process_sample_pair(sample_name, csv_path, fasta_path)
            
            # Step 6: Write sample-specific files sorted by length
            if modified_records:
                sample_output_records = sorted(
                    modified_records, 
                    key=lambda rec: len(rec.seq), 
                    reverse=True
                )
                output_fasta_path = output_dir / f"{sample_name}_annotated.fasta"
                SeqIO.write(sample_output_records, output_fasta_path, "fasta")
                print(f"  -> Successfully wrote sorted FASTA file to {output_fasta_path}")
            
            # Aggregate all records by Taxon
            for rec in modified_records:
                # Use the clean_taxon attribute added in process_sample_pair
                taxon_records_map[rec.clean_taxon].append(rec)
            
        else:
            print(f"  [Warning] Corresponding FASTA file not found: {fasta_filename}. Skipping sample.")


    # Step 7 (Revised): Write records grouped by Taxon
    print("-" * 40)
    print("Writing grouped FASTA files to BY_TAXON directory...")
    
    total_taxa = 0
    total_records = 0
    
    for taxon_name, records in taxon_records_map.items():
        total_taxa += 1
        total_records += len(records)
        
        # Sort records for the current Taxon by length (bigger to smaller)
        records.sort(key=lambda rec: len(rec.seq), reverse=True)
        
        # Write to the file named after the Taxon
        taxon_fasta_path = taxon_output_dir / f"{taxon_name}.fasta"
        SeqIO.write(records, taxon_fasta_path, "fasta")
        
        print(f"  -> Wrote {len(records)} records for Taxon: {taxon_name}.fasta")
        
    print("-" * 40)
    print(f"Processing complete. Wrote files for {total_taxa} unique Taxa, containing {total_records} records.")


if __name__ == '__main__':
    main()

# Example usage:
# python fasta_processor.py /path/to/INPUT_DIR /path/to/OUTPUT_DIR
