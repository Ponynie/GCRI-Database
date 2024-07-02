from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.rdmolops import GetAdjacencyMatrix
import pandas as pd
import numpy as np

def inchi_to_adjmatrix(inchi_string):
    mol = Chem.MolFromInchi(inchi_string)
    mol = Chem.AddHs(mol)  
    
    if mol is None:
        print("Failed to create molecule from InChI")
        return
    
    formula = rdMolDescriptors.CalcMolFormula(mol)
    print(f"Molecular formula: {formula}")
    
    adj = GetAdjacencyMatrix(mol)
    
    element = [atom.GetSymbol() for atom in mol.GetAtoms()]
    print(pd.DataFrame(adj, index=element, columns=element))
    
    return np.float32(adj)

inchi_to_adjmatrix("InChI=1S/H2O/h1H2")
