import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops, Draw
import matplotlib.pyplot as plt
import networkx as nx

def get_molecule_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0, 0, 0  # Handle invalid SMILES
    
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    num_rings = mol.GetRingInfo().NumRings()
    
    return num_atoms, num_bonds, num_rings

def calculate_diameter(adj_matrix):
    G = nx.Graph(adj_matrix)
    if nx.is_connected(G):
        return nx.diameter(G)
    else:
        # Calculate diameter of the largest connected component
        largest_cc = max(nx.connected_components(G), key=len)
        subgraph = G.subgraph(largest_cc)
        return nx.diameter(subgraph)

# Load the dataset
file_name = "data/NP-KV-ALL-G-C.csv"
data = pd.read_csv(file_name)

# Calculate the number of atoms, bonds, and rings for each molecule
data['num_atoms'], data['num_bonds'], data['num_rings'] = zip(*data['smiles'].apply(get_molecule_properties))

# Calculate the sum of atoms and bonds for each molecule
data['atom_bond_sum'] = data['num_atoms'] + data['num_bonds']

# Identify the molecule with the largest atom-bond sum
largest_molecule = data.loc[data['atom_bond_sum'].idxmax()]
largest_smiles = largest_molecule['smiles']

# Visualize the largest molecule
largest_mol = Chem.MolFromSmiles(largest_smiles)
img = Draw.MolToImage(largest_mol)
plt.imshow(img)
plt.title('Largest Molecule by Atom and Bond Sum')
plt.axis('off')
plt.show()

# Calculate the diameter for each molecule's adjacency matrix
data['diameter'] = data['smiles'].apply(lambda s: calculate_diameter(rdmolops.GetAdjacencyMatrix(Chem.MolFromSmiles(s))) if Chem.MolFromSmiles(s) else 0)

# Identify the molecule with the largest diameter
largest_diameter_molecule = data.loc[data['diameter'].idxmax()]
largest_diameter_smiles = largest_diameter_molecule['smiles']

# Visualize the molecule with the largest diameter
largest_diameter_mol = Chem.MolFromSmiles(largest_diameter_smiles)
img_diameter = Draw.MolToImage(largest_diameter_mol)
plt.imshow(img_diameter)
plt.title('Molecule with the Largest Diameter')
plt.axis('off')
plt.show()

# Display statistics
print("Statistics of Number of Atoms:")
print(data['num_atoms'].describe())

print("\nStatistics of Number of Bonds:")
print(data['num_bonds'].describe())

print("\nStatistics of Number of Rings:")
print(data['num_rings'].describe())

print(f"\nLargest Molecule by Atom and Bond Sum:")
print(f"  SMILES: {largest_smiles}")
print(f"  Number of Atoms: {largest_molecule['num_atoms']}")
print(f"  Number of Bonds: {largest_molecule['num_bonds']}")
print(f"  Number of Rings: {largest_molecule['num_rings']}")
print(f"  Atom-Bond Sum: {largest_molecule['atom_bond_sum']}")

print(f"\nMolecule with the Largest Diameter:")
print(f"  SMILES: {largest_diameter_smiles}")
print(f"  Number of Atoms: {largest_diameter_molecule['num_atoms']}")
print(f"  Diameter: {largest_diameter_molecule['diameter']}")