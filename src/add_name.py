import os
import pandas as pd

# Paths
large_file_path = "data/data.csv"  # The large CSV file
small_files_folder = "data/agg_by_Inchi"  # Folder containing smaller CSV files

# Step 1: Load the large CSV file and create a dictionary for InChI to Name mapping
print("Loading the large reference file...")
large_df = pd.read_csv(large_file_path)

# Ensure the required columns exist
if "inChI" not in large_df.columns or "Name" not in large_df.columns:
    raise ValueError("Large file must contain 'inChI' and 'Name' columns.")

# Create a mapping dictionary (drop duplicates to avoid conflicts)
# Ensure we keep the first occurrence of each inChI key
inchi_to_name = large_df.drop_duplicates(subset="inChI", keep="first").set_index("inChI")["Name"].to_dict()

# Step 2: Process each smaller CSV file
if not os.path.exists(small_files_folder):
    raise FileNotFoundError(f"Folder '{small_files_folder}' does not exist.")

for file_name in os.listdir(small_files_folder):
    if file_name.endswith(".csv"):
        file_path = os.path.join(small_files_folder, file_name)
        print(f"Processing file: {file_path}")

        # Read the smaller CSV file
        small_df = pd.read_csv(file_path)

        # Ensure the required column exists
        if "inChI" not in small_df.columns:
            print(f"Skipping {file_name}: Missing 'inChI' column.")
            continue

        # Step 3: Add the Name column by mapping from the dictionary
        small_df["Name"] = small_df["inChI"].map(inchi_to_name)

        # Step 4: Save the updated file
        small_df.to_csv(file_path, index=False)
        print(f"Updated file saved: {file_path}")