import os
import pandas as pd
from sklearn.model_selection import train_test_split

# Base directories
data_dir = 'data'
splitted_dir = 'data/splitted'

# Function to process and split CSV files
def process_and_split_csv(file_path, output_dir, suffix='-TEST'):
    # Load the CSV file
    data = pd.read_csv(file_path)

    # Perform stratified train-test split
    train, test = train_test_split(
        data,
        test_size=0.1,
        stratify=data['stratify_group'],
        random_state=42
    )

    # Determine output file names
    base_name = os.path.basename(file_path)
    train_output_path = os.path.join(output_dir, base_name)
    test_output_path = os.path.join(output_dir, base_name.replace('.csv', f'{suffix}.csv'))

    # Save the split datasets
    os.makedirs(output_dir, exist_ok=True)
    train.to_csv(train_output_path, index=False)
    test.to_csv(test_output_path, index=False)

    print(f"Processed {file_path}. Train set: {train_output_path}, Test set: {test_output_path}")

# Process each specific file
for file_name in os.listdir(data_dir):
    if file_name.endswith('.csv') and '-G-C-P' in file_name:  # Match the required naming pattern
        prefix = '-'.join(file_name.split('-')[:3])  # Extract the prefix, e.g., "NP-KV-ISO"
        file_path = os.path.join(data_dir, file_name)
        split_output_dir = os.path.join(splitted_dir, prefix)  # Use mapped directory
        process_and_split_csv(file_path, split_output_dir)

print("Processing complete.")