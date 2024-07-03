import os
import requests

def download_file(url, dest_directory):
    # Ensure the destination directory exists
    if not os.path.exists(dest_directory):
        os.makedirs(dest_directory)

    # Extract the filename from the URL
    filename = os.path.basename(url)

    # Create the full file path
    file_path = os.path.join(dest_directory, filename)

    # Download the file
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"File downloaded successfully and saved to {file_path}")
    else:
	print(f"Failed to download file. Status code: {response.status_code}")

# URL of the file to download
file_url = 'https://humantfs.ccbr.utoronto.ca/download/v_1.01/DBD_Alignments_v_1.01.zip'

# Destination directory
destination_directory = '/home/bg171/Project/dbds'

# Call the download function
download_file(file_url, destination_directory)
