import os

# Function to extract TAZ2 domain from EP300.domains.txt
def extract_taz2_domain(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    taz2_line = None
    for line in lines:
        if '1660-1840' in line:
            taz2_line = '1660-1840'
            break
    return taz2_line

# Function to process activation domain files
def process_activation_domain_files(taz2_domain, directory):
    processed_lines = []
    for filename in os.listdir(directory):
        if filename.endswith('.fa'):
            with open(os.path.join(directory, filename), 'r') as file:
                lines = file.readlines()
            tf_code = lines[0].strip().lstrip('>')
            ad_locations = []
            for line in lines[1:]:
                if 'Location' in line:
                    ad_locations.append(line.strip().split(": ")[1])
            if ad_locations:
                ad_locations_str = ','.join(ad_locations)
                processed_lines.append(f"{tf_code},{taz2_domain};{tf_code},{ad_locations_str}")
    return '\n'.join(processed_lines)

# Define the path to the EP300.domains.txt file and the directory with .fa files
ep300_domains_file_path = '/home/bg171/Project/p300/EP300.domains.txt'
activation_domain_directory = '/home/bg171/Project/ADs/AD_locations'

# Extract the TAZ2 domain
taz2_domain = extract_taz2_domain(ep300_domains_file_path)

# Process the activation domain files
processed_output = process_activation_domain_files(taz2_domain, activation_domain_directory)

# Write the processed output to a new file
with open('processed_activation_domains.txt', 'w') as file:
    file.write(processed_output)

processed_output
