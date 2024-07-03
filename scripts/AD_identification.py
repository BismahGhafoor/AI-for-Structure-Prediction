import pandas as pd
import os
from Bio import SeqIO
import openpyxl

# Function to find activation domains in a TF sequence
def find_activation_domains(tf_seq, tad_df):
    matches = []
    tf_seq_str = str(tf_seq)
    
    for idx, row in tad_df.iterrows():
        gene = row['Gene']
        fragment = row['Fragment']
        tad_seq = row['Sequence']

        start = tf_seq_str.find(tad_seq)
        if start != -1:
            end = start + len(tad_seq)
            matches.append((gene, fragment, start, end, tad_seq))
    
    return matches

# Read the tAD-seq sheet
file_path = '/home/bg171/Project/ADs/mmc2.xlsx'
tad_seq_data = pd.read_excel(file_path, sheet_name='tAD-seq', engine='openpyxl')

# Directory containing the TF sequences
input_dir = '/home/bg171/Project/FastaFiles'
output_dir = '/home/bg171/Project/ADs/AD_locations'
no_activation_domain_file = os.path.join(output_dir, 'no_activation_domain.txt')

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# List to store TFs with no activation domain
no_activation_domain_list = []

# Process each TF sequence file
for file_name in os.listdir(input_dir):
    if file_name.endswith('.fa'):
        file_path = os.path.join(input_dir, file_name)
        with open(file_path, 'r') as file:
            header = file.readline().strip()
            tf_seq = file.read().replace("\n", "")

        matches = find_activation_domains(tf_seq, tad_seq_data)

        output_file_path = os.path.join(output_dir, file_name)
        with open(output_file_path, 'w') as output_file:
            output_file.write(header + "\n")
            if matches:
                for match in matches:
                    gene, fragment, start, end, tad_seq = match
                    output_file.write(f"Location: {start}-{end}\n")
            else:
                output_file.write("No activation domain detected\n")
                no_activation_domain_list.append(file_name)

# Write the no activation domain list to a file
with open(no_activation_domain_file, 'w') as file:
    for tf_file in no_activation_domain_list:
        file.write(tf_file + "\n")

print("Processing complete. Check the output directory for results.")
