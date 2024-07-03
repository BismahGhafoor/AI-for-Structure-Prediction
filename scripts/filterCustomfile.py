import re

# Define a function to filter ranges
def filter_ranges(line):
    # Split the line into parts
    parts = line.strip().split(';')
    filtered_parts = []

    for part in parts:
        # Find the first comma, which separates the protein from the ranges
        first_comma_index = part.find(',')
        protein = part[:first_comma_index]
        ranges = part[first_comma_index + 1:]

        # Split ranges by comma and filter out invalid ranges
        range_list = ranges.split(',')
        filtered_ranges = [r for r in range_list if not re.match(r'1-(\d{1,2}|100)$', r)]

        if filtered_ranges:
            filtered_parts.append(f"{protein},{','.join(filtered_ranges)}")

    return ';'.join(filtered_parts)

# Read the input file
with open('custom.txt', 'r') as file:
    lines = file.readlines()

# Process each line
filtered_lines = [filter_ranges(line) for line in lines]

# Remove lines where nothing is left after filtering
filtered_lines = [line for line in filtered_lines if ';' in line and line.strip().split(';')[1]]

# Write the output to a new file or overwrite the existing file
with open('custom_filtered.txt', 'w') as file:
    file.write('\n'.join(filtered_lines))

print("Processing complete. Check the custom_filtered.txt file for results.")
