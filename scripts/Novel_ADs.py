import os
import requests
from Bio import SeqIO
import time
import csv
import pickle
from io import StringIO
import af2plots
print("af2plots module imported successfully!")
import sys
import pickle
import io
import os

print(sys.executable)

class DummyJaxArray:
    def __init__(self, *args, **kwargs):
        self.shape = args[0] if args else None
        self.dtype = kwargs.get('dtype', None)

    def __getattr__(self, name):
        return lambda *args, **kwargs: self

    def __array__(self):
        import numpy as np
        if self.shape:
            return np.zeros(self.shape)
        return np.array([])

class CustomUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "jax._src.device_array":
            return DummyJaxArray
        return super().find_class(module, name)

def custom_load_pickle(file):
    if isinstance(file, (str, bytes, os.PathLike)):
        try:
            with open(file, 'rb') as f:
                return CustomUnpickler(f).load()
        except pickle.UnpicklingError as e:
            print(f"Error unpickling file {file}: {e}")
            # Handle or re-raise the error as appropriate
            return None
    elif isinstance(file, io.IOBase):
        try:
            return CustomUnpickler(file).load()
        except pickle.UnpicklingError as e:
            print(f"Error unpickling file {file.name if hasattr(file, 'name') else '<stream>'}: {e}")
            return None
    else:
        raise TypeError(f"Unsupported file type: {type(file)}")


# Monkey-patch the pickle.load function in af2plots
import af2plots.plotter
af2plots.plotter.pickle.load = custom_load_pickle

# Paths
fasta_dir = "/scratch/alice/b/bg171/FinalProject/FastaFiles"
dbd_dir = "/scratch/alice/b/bg171/FinalProject/dbd_locations"
pkl_base_dir = "/scratch/alice/b/bg171/FinalProject/Predictions"  # Base directory containing directories with pkl files
output_file = "Novel_ADs.csv"

def download_alphamissense_data(uniprot_id):
    print(f"Downloading AlphaMissense data for {uniprot_id}")
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-aa-substitutions.csv"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        csv_content = StringIO(response.text)
        reader = csv.DictReader(csv_content)

        alphamissense_data = {}
        for row in reader:
            variant = row['protein_variant']
            position = int(variant[1:-1])  # Extract position from variant (e.g., "M1A" -> 1)
            am_pathogenicity = float(row['am_pathogenicity'])
            alphamissense_data[position] = am_pathogenicity

        print(f"Successfully downloaded data for {uniprot_id}")
        return alphamissense_data
    except requests.exceptions.RequestException as e:
        print(f"Error downloading data for {uniprot_id}: {e}")
    return None

def parse_location_file(file_path):
    locations = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Location:"):
                start, end = map(int, line.split(":")[1].strip().split("-"))
                locations.append((start, end))
    print(f"Parsed locations from {file_path}: {locations}")
    return locations

def load_plddt_data_from_pkl(pkl_file_path):
    """
    Load pLDDT scores from the .pkl file.
    """
    with open(pkl_file_path, 'rb') as f:
        data = custom_load_pickle(f)
        plddt_scores = {i + 1: plddt for i, plddt in enumerate(data['plddt'])}
    return plddt_scores

def find_pkl_file_for_uniprot(pkl_base_dir, uniprot_id):
    """
    Find the correct .pkl file for a given Uniprot ID within the pkl_base_dir.
    """
    # Traverse the base directory
    for root, dirs, files in os.walk(pkl_base_dir):
        # Check each subdirectory
        for dir_name in dirs:
            # Check if the subdirectory name contains the Uniprot ID
            if uniprot_id in dir_name:
                # Construct the full path to the subdirectory
                subdir_path = os.path.join(root, dir_name)
                # Construct the expected path to the .pkl file
                pkl_file_path = os.path.join(subdir_path, "result_model_1_multimer_v3_pred_0.pkl")
                # Check if the .pkl file exists
                if os.path.exists(pkl_file_path):
                    print(f"Found .pkl file for {uniprot_id}: {pkl_file_path}")
                    return pkl_file_path
    # If not found, log a warning and return None
    print(f"No .pkl file found for {uniprot_id}")
    return None
def identify_regions_with_plddt(alphamissense_data, plddt_data, threshold_high=0.564, threshold_low=0.2, plddt_high=60, plddt_low=50):
    """
    Identify regions with high AFmissense scores flanked by low AFmissense scores
    and check for the required pLDDT score patterns.
    """
    print("Identifying regions based on AFmissense and pLDDT score criteria")
    regions = []
    keys = sorted(alphamissense_data.keys())

    print(f"AFmissense data keys: {keys[:10]}")  # Print first 10 keys for reference
    print(f"First 10 pLDDT values: {[plddt_data[k] for k in keys[:10]]}")

    in_high_region = False
    start = None
    end = None
    
    for i in range(1, len(keys) - 1):
        prev_score = alphamissense_data[keys[i - 1]]
        current_score = alphamissense_data[keys[i]]
        next_score = alphamissense_data[keys[i + 1]]

        prev_plddt = plddt_data.get(keys[i - 1], 0)
        current_plddt = plddt_data.get(keys[i], 0)
        next_plddt = plddt_data.get(keys[i + 1], 0)

        print(f"Position {keys[i]}: AFmissense={current_score}, pLDDT={current_plddt}")

        # Adjusted condition: focus on identifying high AFmissense regions, possibly flanked by stable pLDDT
        if (current_score >= threshold_high and 
            (prev_score <= threshold_low or next_score <= threshold_low) and
            (prev_plddt >= plddt_high or next_plddt >= plddt_high)):
            if not in_high_region:
                start = keys[i]
                in_high_region = True
        elif in_high_region and (current_score < threshold_high or next_score <= threshold_low):
            end = keys[i]
            if start is not None and end is not None:
                avg_afmissense = sum(alphamissense_data[j] for j in range(start, end + 1)) / (end - start + 1)
                avg_plddt = sum(plddt_data[j] for j in range(start, end + 1)) / (end - start + 1)
                regions.append((start, end, avg_afmissense, avg_plddt))
            in_high_region = False
            start = None
            end = None

    print(f"Identified regions: {regions}")
    return regions

def main():
    with open(output_file, 'w') as out_f:
        out_f.write("TF\tUniprotID\tRegion\tAverage_AFmissense\tAverage_pLDDT\n")

        for filename in os.listdir(fasta_dir):
            if filename.endswith(".fa"):
                tf_name, uniprot_id = filename.split("_")
                uniprot_id = uniprot_id[:-3]  # Remove .fa
                print(f"\nProcessing {filename}")

                # Read FASTA sequence
                fasta_path = os.path.join(fasta_dir, filename)
                try:
                    with open(fasta_path, 'r') as fasta_file:
                        content = fasta_file.read()
                        print(f"Content of {fasta_path}:")
                        print(content)
                        fasta_file.seek(0)  # Reset file pointer to the beginning
                        seq_record = next(SeqIO.parse(fasta_file, "fasta"))
                        print(f"Parsed sequence: {seq_record.seq[:50]}...")  # Print first 50 characters
                except StopIteration:
                    print(f"Error: Empty or invalid FASTA file: {fasta_path}")
                    continue

                # Get DBD locations
                dbd_file = os.path.join(dbd_dir, f"{filename}.dbd_locations.txt")
                if not os.path.exists(dbd_file):
                    print(f"Warning: DBD file not found: {dbd_file}")
                    dbd_locations = []
                else:
                    dbd_locations = parse_location_file(dbd_file)

                # Find the corresponding pkl file
                pkl_file = find_pkl_file_for_uniprot(pkl_base_dir, uniprot_id)
                if not pkl_file:
                    print(f"Warning: PKL file not found for Uniprot ID {uniprot_id}")
                    continue

                # Load pLDDT scores from the pkl file
                plddt_data = load_plddt_data_from_pkl(pkl_file)

                # Download AlphaMissense data
                alphamissense_data = download_alphamissense_data(uniprot_id)

                if alphamissense_data:
                    # Identify regions based on AFmissense and pLDDT criteria
                    regions = identify_regions_with_plddt(alphamissense_data, plddt_data)

                    for start, end, avg_afmissense, avg_plddt in regions:
                        out_f.write(f"{tf_name}\t{uniprot_id}\t{start}-{end}\t{avg_afmissense:.4f}\t{avg_plddt:.4f}\n")
                else:
                    print(f"Warning: No AlphaMissense data found for {uniprot_id}")

                # Add a delay to avoid overwhelming the API
                time.sleep(1)

if __name__ == "__main__":
    main()
