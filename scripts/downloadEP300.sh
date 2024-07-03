#!/bin/bash

output_dir="/home/bg171/Project/p300" # Define the directory where you want to save the FASTA files
mkdir -p "$output_dir"  # Create the directory if it doesn't exist

log_file="${output_dir}/NoUniProtID_log.txt" # Log file to record entries with NoUniProtID
touch "$log_file" # Create the log file if it doesn't exist

base_url="https://rest.uniprot.org/uniprotkb/search?format=fasta&query=" # Base URL for the UniProt API, specifying the FASTA format for proteins of Homo sapiens

gene_name="EP300" # Define the gene name

query_url="${base_url}(gene:${gene_name}+AND+reviewed:true+AND+organism_id:9606)" # Construct the full query URL for the current gene, ensuring reviewed entries
temp_output_file="${output_dir}/${gene_name}_temp.fasta" # Define a temporary output file path for the current gene's FASTA
echo "Downloading FASTA for ${gene_name}..."
curl -s "${query_url}" -o "${temp_output_file}" # Use curl to fetch the FASTA data and save it to the temporary output file

# Check if the file has content, if not report as no data found
if ! [ -s "$temp_output_file" ]; then
    echo "No data found for ${gene_name}, removing empty file."
    rm "$temp_output_file"
else
    # Process each FASTA record separately
    awk -v output_dir="$output_dir" -v gene_name="$gene_name" -v log_file="$log_file" '
    BEGIN { RS=">"; FS="\n" }
    NR > 1 {
	header = ">" $1
        seq = ""
        for (i=2; i<=NF; i++) {
            seq = seq $i
        }
	if (header ~ /^>sp\|/) {
            split(header, arr, "|")
            uniprot_id = arr[2]
        } else {
            uniprot_id = "NoUniProtID"
            print gene_name >> log_file
        }
	final_output_file = output_dir "/" gene_name "_" uniprot_id ".fa"
        if (uniprot_id == "NoUniProtID") {
            print header > final_output_file
        } else {
            print ">" uniprot_id > final_output_file
        }
	print seq >> final_output_file
    }
    ' "$temp_output_file"

    # Remove the temporary file
    rm "$temp_output_file"
fi

echo "Download complete." # Print 'Download complete' when all FASTA files have been successfully downloaded
