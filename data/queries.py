import pandas as pd
import os
from rdkit import Chem

ri_types = ['Van Den Dool and Kratz', 'Normal alkane', "Kovats'", "Lee's"]
phase_types = ['non-polar', 'polar']
temperature_modes = ['temperature ramp', 'custom temperature program', 'isothermal']

def ensure_relative_path(path):
    current_dir = os.path.dirname(__file__)
    return os.path.join(current_dir, '..', path)

def separate_detail_column():
    
    combined_data = pd.read_csv('data/data-combined.csv')

    # Count the occurrences of unique values in the "Details" column
    details_counts = combined_data["Details"].value_counts()

    # Print the counts
    print("Unique values in the 'Details' column and their counts:")
    for detail, count in details_counts.items():
        print(f"{detail}: {count}")

    df = combined_data

    # Define regular expressions to extract RI Type, Column Type, and Temperature Mode
    ri_type_pattern = r"([\w'\s]+) RI"
    column_type_pattern = r", (non-polar|polar) column"
    temperature_mode_pattern = r", (isothermal|temperature ramp|custom temperature program)"

    # Extract RI Type, Column Type, and Temperature Mode from the "Details" column
    df["RI Type"] = df["Details"].str.extract(ri_type_pattern, expand=False)
    df["Phase Polarity"] = df["Details"].str.extract(column_type_pattern, expand=False)
    df["Temperature Mode"] = df["Details"].str.extract(temperature_mode_pattern, expand=False)

    # Drop the "Details" column
    df = df.drop("Details", axis=1)

    # Save the modified DataFrame to a new CSV file
    df.to_csv('data/data.csv', index=False)

def queries(**kwargs):
    
    # Read the CSV data
    data_path = ensure_relative_path('data/data.csv')
    df = pd.read_csv(data_path)

    # Create a dictionary mapping keyword arguments to column names
    column_mapping = {
        'RI_Type': 'RI Type',
        'Phase_Polarity': 'Phase Polarity',
        'Temperature_Mode': 'Temperature Mode',
        'Active_Phase': 'Active phase',
        'Column_Type': 'Column type',
    }

    # Filter the DataFrame based on the specified chromatographic conditions
    filter_conditions = pd.Series([True] * len(df), index=df.index)
    for arg_name, value in kwargs.items():
        column_name = column_mapping.get(arg_name, arg_name)
        if column_name in df.columns:
            if '|' in str(value):
                filter_conditions &= df[column_name].isin(value.split('|'))
            else:
                filter_conditions &= (df[column_name] == value)

    filtered_df = df[filter_conditions]

    # Define the constant group_by_columns (unique identifiers for each compound)
    group_by_columns = ['inChI']

    # Aggregate the duplicates based on retention index (I) column
    aggregate_column = 'I'
    aggregated_df = filtered_df.groupby(group_by_columns)[aggregate_column].mean().reset_index()

    return aggregated_df

def view_data_number():
    s = 0
    for ri_type in ri_types:
        for phase_type in phase_types:
            for temperature_mode in temperature_modes:
                print("-----------------------------------------------------------")
                print(f"RI Type: {ri_type}, Phase Polarity: {phase_type}, Temperature Mode: {temperature_mode}")
                print(len(queries(RI_Type=ri_type, Phase_Polarity=phase_type, Temperature_Mode=temperature_mode)))
                s += len(queries(RI_Type=ri_type, Phase_Polarity=phase_type, Temperature_Mode=temperature_mode))
                print("-----------------------------------------------------------")
    print(s)

def prepare_traintest_data():
    df = queries(RI_Type="Van Den Dool and Kratz|Normal alkane|Kovats'", Phase_Polarity="non-polar")
    #df.drop(['molecularFormula', 'inChIKey'], axis=1, inplace=True)

    def inchi_to_smiles(inchi):
        mol = Chem.MolFromInchi(inchi)
        if mol:
            return Chem.MolToSmiles(mol)
        else:
            return None

    df['SMILES'] = df['inChI'].apply(inchi_to_smiles)
    
    print('Number of queries compounds with SMILES:', df['SMILES'].notnull().sum())

    df.drop('inChI', axis=1, inplace=True)
    df.rename(columns={'I': 'ri'}, inplace=True)
    df.rename(columns={'SMILES': 'smiles'}, inplace=True)
    df.dropna(subset=['smiles'], inplace=True)
    df.to_csv(ensure_relative_path('data/NonPolarRI.csv'), index=False)
    
def uniqueness_test():
    df = pd.read_csv(ensure_relative_path('data/NonPolarRI.csv'))
    print('Number of unique compounds:', len(df['smiles'].unique()))
    print('Number of compounds:', len(df))
    
if __name__ == '__main__':
    uniqueness_test()