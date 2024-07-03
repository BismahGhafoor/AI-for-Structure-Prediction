#!/bin/bash

output_dir="/home/bg171/Project/fasta_files" #  Directory to save the FASTA files
mkdir -p "$output_dir"  # Create the directory if it doesn't exist

base_url="https://rest.uniprot.org/uniprotkb/search?format=fasta&query=" # Base URL for the UniProt API, specifying the FASTA format for proteins of Homo sapiens

# Loop through each gene name in the genenames.txt file
while IFS= read -r gene_name; do
    query_url="${base_url}(gene:${gene_name}+AND+reviewed:true+AND+organism_id:9606)" # Construct the full query URL for the current gene, ensuring reviewed entries
    output_file="${output_dir}/${gene_name}.fasta" # Output file path for the current gene's FASTA
    echo "Downloading FASTA for ${gene_name}..."
    curl -s "${query_url}" -o "${output_file}" # Use curl to fetch the FASTA data and save it to the output file

    # Check if the file has content, if not report as no data found
    if ! [ -s "$output_file" ]; then
        echo "No data found for ${gene_name}, removing empty file."
        rm "$output_file"
    fi
done < "genenames.txt"

echo "Download complete." # Print 'Download complete' when all FASTA files have been successfully downloaded
