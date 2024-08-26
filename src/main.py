from rdkit import Chem

water = Chem.MolFromSmiles('[OH2]')
a = set(atom.GetSymbol() for atom in water.GetAtoms())
print(a)