import json

# Function to read the Rochester correction text files
def read_rochester_file(file_path):
    corrections = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue  # Skip comments and empty lines
            parts = line.split()
            corrections.append([float(p) for p in parts])
    return corrections

# Function to transform the corrections into a JSON format
def transform_to_json(corrections):
    json_data = {
        "version": 1,
        "name": "rochester_correction",
        "description": "Rochester correction for muons",
        "inputs": [
            {"name": "pt", "type": "real"},
            {"name": "eta", "type": "real"},
            {"name": "phi", "type": "real"},
            {"name": "charge", "type": "real"},
        ],
        "output": {"name": "correction", "type": "real"},
        "data": corrections
    }
    return json_data

# Read the corrections from the text file
file_path = 'RoccoR2018.txt'
corrections = read_rochester_file(file_path)

# Transform to JSON
json_data = transform_to_json(corrections)

# Save the JSON data to a file
with open('rochester_correction.json', 'w') as json_file:
    json.dump(json_data, json_file, indent=4)

print("Rochester correction data has been converted to JSON and saved.")