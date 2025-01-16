import pandas as pd

# File names and corresponding condition labels
files_conditions = {
    "data/NP-KV-ISO-G-C-P.csv": "NP-KV",
    "data/P-KV-ISO-G-C-P.csv": "P-KV",
    "data/NP-VAN-RAMP-G-C-P.csv": "NP-VAN",
    "data/P-VAN-RAMP-G-C-P.csv": "P-VAN"
}

# Create an empty dataframe to hold combined data
combined_df = pd.DataFrame()

# Loop through each file and read data
for file_name, condition in files_conditions.items():
    df = pd.read_csv(file_name)  # Read the CSV file
    df = df.rename(columns={'ri': f'{condition}_value'})  # Rename 'ri' to a condition-specific column name
    df[condition] = 1  # Binary column indicating the presence of the RI value (always 1 if present)

    # Drop 'group' and 'stratify_group' columns as they are not needed for now
    df = df.drop(columns=['group', 'stratify_group'], errors='ignore')

    if combined_df.empty:
        combined_df = df[['smiles', f'{condition}_value', condition]]
    else:
        combined_df = pd.merge(combined_df, df[['smiles', f'{condition}_value', condition]], on='smiles', how='outer')  # This merges based on 'smiles', adding all unique molecules with NaN-filled columns where necessary for missing RI values.

# Fill missing binary indicators and RI values
for condition in files_conditions.values():
    combined_df[condition] = combined_df[condition].fillna(0)  # Fill binary indicator with 0 if molecule not in condition using assignment to avoid inplace warning
    combined_df[f'{condition}_value'] = combined_df[f'{condition}_value'].fillna(pd.NA)  # Fill missing RI values with NaN (pd.NA) instead of 0

# Save the combined dataframe to CSV
combined_df.to_csv("data/ALL-KVVAN-ISORAMP.csv", index=False)

print("File saved as ALL-KVVAN-ISORAMP.csv")