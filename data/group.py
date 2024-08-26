import pandas as pd
from rdkit import Chem
from rdkit.Chem import Fragments
import inspect
from collections import Counter

def get_fragment_functions():
    return {name: func for name, func in inspect.getmembers(Fragments) if name.startswith('fr_')}

def count_fragments(mol, fragment_funcs: dict):
    return {name: func(mol) for name, func in fragment_funcs.items()}

def assign_group(fragment_counts):
    max_count = max(fragment_counts.values())
    if max_count == 0:
        return 'unknown'
    return max(fragment_counts, key=fragment_counts.get)

def add_molecular_group(df):
    fragment_funcs = get_fragment_functions()
    
    def process_row(row):
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is None:
            return pd.Series({'group': 'invalid_smiles'})
        fragment_counts = count_fragments(mol, fragment_funcs)
        group = assign_group(fragment_counts).replace('fr_', '')
        return pd.Series({'group': group})
    
    result = df.apply(process_row, axis=1)
    return pd.concat([df, result], axis=1)

def handle_rare_groups(df, min_samples=2):
    df['stratify_group'] = df['group']
    
    group_counts = Counter(df['group'])
    rare_groups = {group for group, count in group_counts.items() if count < min_samples}
    
    if len(rare_groups) == 1:
        # If there's only one rare instance, assign it to the second least common group
        second_least_common = sorted(group_counts.items(), key=lambda x: x[1])[1][0]
        df.loc[df['group'].isin(rare_groups), 'stratify_group'] = second_least_common
    else:
        df['stratify_group'] = df['group'].apply(lambda x: 'rare' if x in rare_groups else x)
    
    return df

# Load your data
df = pd.read_csv('data/NP-LRI-RAMP-C.csv')

# Add molecular groups
df_with_groups = add_molecular_group(df)

# Remove invalid SMILES
df_with_groups = df_with_groups[df_with_groups['group'] != 'invalid_smiles']

# Handle rare groups
df_with_groups = handle_rare_groups(df_with_groups, min_samples=10)

# Save the result with groups
df_with_groups.to_csv('data/NP-LRI-RAMP-G-C.csv', index=False)
