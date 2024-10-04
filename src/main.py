import pandas as pd

# Define the file path and sheet names
file_path = 'data/data-GCxGC.xlsx'  # Update this to the correct path if necessary
sheets = ["Training set", "Test set", "External validation set"]
output_files = ["data/data-GCxGC_train.csv", 
                "data/data-GCxGC_test.csv", 
                "data/data-GCxGC_validation.csv"]

# Define the columns to keep and the new column names
columns_to_keep = ["SMILES", "Alkane RI", "PEG-2I (values for compounds with -999 could not be calculated)"]
new_column_names = {"SMILES": "smiles", "Alkane RI": "ri", 
                    "PEG-2I (values for compounds with -999 could not be calculated)": "ri_2"}

# Load each sheet, filter the columns, rename them, remove rows where PEG-2I equals -999, then save to CSV
for sheet, output_file in zip(sheets, output_files):
    # Load the sheet
    df = pd.read_excel(file_path, sheet_name=sheet)
    
    # Filter the columns
    df_filtered = df[columns_to_keep]
    
    # Rename columns
    df_filtered = df_filtered.rename(columns=new_column_names)
    
    # Remove rows where 'ri_2' equals -999
    df_filtered = df_filtered[df_filtered["ri_2"] != -999]
    
    # Save to CSV
    df_filtered.to_csv(output_file, index=False)

print("CSV files created and columns renamed successfully!")
a = "test"