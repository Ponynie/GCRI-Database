import pandas as pd

# Define the file path and sheet names
file_path = 'data/data-GCxGC.xlsx'  # Update this to the correct path if necessary
sheets = ["Training set", "Test set", "External validation set"]
output_files = ["data/data-GCxGC_train.csv", 
                "data/data-GCxGC_test.csv", 
                "data/data-GCxGC_validation.csv"]

# Define the columns to keep
columns_to_keep = ["SMILES", "Alkane RI", "PEG-2I (values for compounds with -999 could not be calculated)"]

# Load each sheet, filter the columns, and remove rows where PEG-2I equals -999, then save to CSV
for sheet, output_file in zip(sheets, output_files):
    # Load the sheet
    df = pd.read_excel(file_path, sheet_name=sheet)
    
    # Filter the columns
    df_filtered = df[columns_to_keep]
    
    # Remove rows where 'PEG-2I' equals -999
    df_filtered = df_filtered[df_filtered[columns_to_keep[2]] != -999]
    
    # Save to CSV
    df_filtered.to_csv(output_file, index=False)

print("CSV files created successfully!")