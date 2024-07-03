import os

def extract_non_dbd_regions(fa_directory, dbd_directory, output_file, p300_domains_file, valid_dbd_log):
    results = []

    # Read p300 domains
    with open(p300_domains_file, 'r') as p300_file:
        p300_lines = p300_file.readlines()
        p300_id = p300_lines[0].strip().replace('>', '')
        p300_locations = p300_lines[1].split('Location:')[1].strip().replace(' ', '')

    # Read valid DBD log
    with open(valid_dbd_log, 'r') as log_file:
        valid_files = {line.strip() for line in log_file}

    # List all .fa files in the given directory
    for filename in os.listdir(fa_directory):
        if filename.endswith('.fa') and filename in valid_files:
            fa_filepath = os.path.join(fa_directory, filename)
            dbd_filepath = os.path.join(dbd_directory, filename + '.dbd_locations.txt')

            # Extract the sequence length from the .fa file
            with open(fa_filepath, 'r') as fa_file:
                lines = fa_file.readlines()
                uniprot_id = lines[0].strip().replace('>', '')
                sequence = ''.join(lines[1:]).replace('\n', '')
                sequence_length = len(sequence)

            # Extract the DBD locations
            non_dbd_regions = []
            dbd_ranges = []
            if os.path.isfile(dbd_filepath):
                with open(dbd_filepath, 'r') as dbd_file:
                    lines = dbd_file.readlines()
                    for line in lines:
                        if line.startswith('Location:'):
                            locations = line.split('Location:')[1].strip()
                            for location in locations.split(','):
                                try:
                                    dbd_start, dbd_end = map(int, location.split('-'))
                                    dbd_ranges.append((dbd_start, dbd_end))
                                except ValueError:
                                    print(f"Invalid location format in file {dbd_filepath}: {location}")
                                    continue
                dbd_ranges.sort()  # Ensure the ranges are sorted by start position

                # Calculate non-DBD regions
                previous_end = 0
                for dbd_start, dbd_end in dbd_ranges:
                    if dbd_start > previous_end + 1:
                        non_dbd_regions.append(f"{previous_end + 1}-{dbd_start - 1}")
                    previous_end = dbd_end
                if previous_end < sequence_length:
                    non_dbd_regions.append(f"{previous_end + 1}-{sequence_length}")
            else:
                # If no DBD file is found, consider the entire sequence as non-DBD
                non_dbd_regions.append(f"1-{sequence_length}")

            # Combine the results for this file into one line
            if non_dbd_regions:
                non_dbd_regions_str = ",".join(non_dbd_regions).replace(' ', '')
                result_line = f"{p300_id},{p300_locations};{uniprot_id},{non_dbd_regions_str}"
                results.append(result_line)
    
    # Write the results to the output file
    with open(output_file, 'w') as output:
        for result in results:
            output.write(result + '\n')
    
    # Print the total number of entries
    print(f"Total number of entries: {len(results)}")

# Define the directory paths and output file name
fa_directory = '/home/bg171/Project/FastaFiles'  # Replace with the path to your .fa files directory
dbd_directory = '/home/bg171/Project/dbds/dbd_locations'  # Replace with the path to your .dbd_locations files directory
output_file_name = 'custom.txt'
p300_domains_file = '/home/bg171/Project//p300/EP300.domains.txt'  # Replace with the path to your p300 domains file
valid_dbd_log = '/home/bg171/Project/dbds/dbd_locations/valid_dbd_log.txt'  # Path to the valid DBD log file

# Run the extraction
extract_non_dbd_regions(fa_directory, dbd_directory, output_file_name, p300_domains_file, valid_dbd_log)
