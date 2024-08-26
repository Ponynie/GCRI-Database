import pandas as pd
from rdkit import Chem
from collections import defaultdict

# Load your CSV file
df = pd.read_csv('data/NP-LRI-RAMP-G-C.csv')

# Initialize a dictionary to count unique molecules containing each atom type
atom_molecule_count = defaultdict(int)

# First pass: count the number of molecules containing each atom type
for smiles in df['smiles']:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        unique_atoms = set(atom.GetSymbol() for atom in mol.GetAtoms())
        for atom in unique_atoms:
            atom_molecule_count[atom] += 1

# Filter out atoms that appear in fewer than 50 different molecules
atoms_to_keep = {atom for atom, count in atom_molecule_count.items() if count >= 50}
print(f"Atoms to keep: {atoms_to_keep}")

# Second pass: filter the dataframe to keep only molecules with the desired atom types
def molecule_has_valid_atoms(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        unique_atoms = set(atom.GetSymbol() for atom in mol.GetAtoms())
        return all(atom in atoms_to_keep for atom in unique_atoms)
    return False

filtered_df = df[df['smiles'].apply(molecule_has_valid_atoms)]

# Save the filtered data to a new CSV file
filtered_df.to_csv('data/NP-LRI-RAMP-G-C-P.csv', index=False)
