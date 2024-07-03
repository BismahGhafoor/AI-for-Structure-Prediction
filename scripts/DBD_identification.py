import os
import zipfile
from Bio import SeqIO
import requests

# Paths
zip_file_path = "/home/bg171/Project/dbds/DBD_Alignments_v_1.01.zip"
extracted_dir = "/home/bg171/Project/dbds/DBD_Alignments"
alignment_subdir = "2018_AddAlignments"
alignment_dir = os.path.join(extracted_dir, alignment_subdir)
fasta_dir = "/home/bg171/Project/FastaFiles"
output_dir = "/home/bg171/Project/dbds/dbd_locations"

# Create necessary directories
os.makedirs(extracted_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)

# Comprehensive keywords associated with TF DNA-binding domains (all in lowercase)
tf_dbd_keywords = [
    "homeobox", "helix-turn-helix", "zinc finger", "leucine zipper", 
    "winged helix", "forkhead", "ets", "tea", "gata-type", "nr c4-type",
    "hox", "myb", "rel", "pou", "mhox", "dm", "bhlh", "zf-c2h2", "zf-c4",
    "zf-c3h1", "t-box", "sry", "ctf/nf-i", "nfat", "dna-binding", "dna binding",
    "dna-binding", "dna binding", "nucleic acid binding", "sequence-specific dna binding"
]

# Step 1: Extract the specific subdirectory from the ZIP file
print(f"Extracting files from {zip_file_path}...")
with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
    zip_ref.extractall(extracted_dir)
print(f"Extraction complete. Extracted files are located in {extracted_dir}")

# Function to parse alignment files and extract sequences without gaps, including consensus sequence
def parse_alignment_files(alignment_dir):
    dbd_sequences = {}
    for root, dirs, files in os.walk(alignment_dir):
        for filename in files:
            if filename.endswith(".fa"):
                filepath = os.path.join(root, filename)
                family = os.path.splitext(filename)[0]
                print(f"Processing alignment file: {filepath}")
                with open(filepath, 'r') as file:
                    sequences = []
                    current_sequence = ""
                    for line in file:
                        if line.startswith(">"):
                            if current_sequence:
                                sequences.append(current_sequence.replace("-", ""))
                            current_sequence = ""
                        else:
                            current_sequence += line.strip()
                    if current_sequence:
                        sequences.append(current_sequence.replace("-", ""))
                    if sequences:
                        dbd_sequences[family] = sequences
    return dbd_sequences

# Function to identify DBD locations in a given sequence
def find_dbd_locations(sequence, dbd_sequences):
    locations = []
    for family, sequences in dbd_sequences.items():
        for dbd_seq in sequences:
            start = sequence.find(dbd_seq)
            if start != -1:
                end = start + len(dbd_seq)
                locations.append((start + 1, end))  # Convert to 1-based index
    return locations

# Function to fetch and parse UniProt entry
def fetch_uniprot_entry(sequence_id):
    url = f"https://rest.uniprot.org/uniprotkb/{sequence_id}.json"
    headers = {
        "User-Agent": "Mozilla/5.0",
    }
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
	print(f"Failed to fetch UniProt entry for {sequence_id}: {response.status_code}")
        return None

# Function to extract relevant features based on keywords from UniProt entry
def extract_relevant_features(uniprot_entry, keywords):
    features = []
    if uniprot_entry:
        for feature in uniprot_entry.get('features', []):
            feature_type = feature.get('type', '').lower()
            description = feature.get('description', '').lower()
            note = feature.get('note', '').lower()
            start = feature.get('location', {}).get('start', {}).get('value')
            end = feature.get('location', {}).get('end', {}).get('value')

            # Check if any keyword matches the type, description, or note
            if any(keyword in feature_type for keyword in keywords) or \
               any(keyword in description for keyword in keywords) or \
               any(keyword in note for keyword in keywords):
                features.append((start, end))
    return features

# Function to merge overlapping and adjacent ranges
def merge_ranges(ranges):
    if not ranges:
        return []
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged_ranges = [sorted_ranges[0]]

    for current in sorted_ranges[1:]:
        last = merged_ranges[-1]
        if current[0] <= last[1] + 1:
            merged_ranges[-1] = (last[0], max(last[1], current[1]))
        else:
            merged_ranges.append(current)
    return merged_ranges

# Step 2: Parse alignment files to get DBD sequences
dbd_sequences = parse_alignment_files(alignment_dir)

# Step 3: Process each FASTA file and identify DBD locations
missing_dbd_log = os.path.join(output_dir, "missing_dbd_log.txt")
valid_dbd_log = os.path.join(output_dir, "valid_dbd_log.txt")
with open(missing_dbd_log, 'w') as missing_log, open(valid_dbd_log, 'w') as valid_log:
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(".fa"):
            fasta_path = os.path.join(fasta_dir, fasta_file)
            output_file = os.path.join(output_dir, f"{fasta_file}.dbd_locations.txt")

            with open(output_file, 'w') as out_f:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    sequence = str(record.seq)
                    sequence_length = len(sequence)
                    print(f"Processing sequence {record.id} with length {sequence_length}")
                    dbd_locations = find_dbd_locations(sequence, dbd_sequences)

                    out_f.write(f">{record.id}\n")
                    if dbd_locations:
                        # Filter out full sequence range if it appears
                        filtered_locations = [
                            (start, end) for start, end in dbd_locations if not (start == 1 and end == sequence_length)
                        ]
                        if filtered_locations:
                            merged_locations = merge_ranges(filtered_locations)
                            if merged_locations:
                                start, end = merged_locations[0]
                                locations_str = f"{start}-{end}"
                                out_f.write(f"Location: {locations_str}\n")
                                valid_log.write(f"{fasta_file}\n")
                            else:
                                out_f.write("No valid DBD locations found.\n")
                                missing_log.write(f"{fasta_file}\n")
                        else:
                            out_f.write("No valid DBD locations found.\n")
                            missing_log.write(f"{fasta_file}\n")
                    else:
                        # If no locations found, search UniProt API
                        uniprot_entry = fetch_uniprot_entry(record.id)
                        features = extract_relevant_features(uniprot_entry, tf_dbd_keywords)
                        if features:
                            filtered_features = [
                                (start, end) for start, end in features if not (start == 1 and end == sequence_length)
                            ]
                            if filtered_features:
                                merged_features = merge_ranges(filtered_features)
                                if merged_features:
                                    start, end = merged_features[0]
                                    locations_str = f"{start}-{end}"
                                    out_f.write(f"Location: {locations_str}\n")
                                    valid_log.write(f"{fasta_file}\n")
                                else:
                                    out_f.write("No valid features found in UniProt search results.\n")
                                    missing_log.write(f"{fasta_file}\n")
                            else:
                                out_f.write("No valid features found in UniProt search results.\n")
                                missing_log.write(f"{fasta_file}\n")
                        else:
                            out_f.write("No relevant features found in UniProt search results.\n")
                            missing_log.write(f"{fasta_file}\n")

print("DBD location identification complete. Check missing_dbd_log.txt for missing DBDs and valid_dbd_log.txt for valid DBDs.")
