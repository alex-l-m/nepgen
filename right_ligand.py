import pandas as pd
from rdkit import Chem

# Currently copied from UYOKUR
right_ligands = ['Fc1ccc2c(c1)[Ir]1n3nc(C(F)(F)F)cc3-c3cccc(n3->1)O2']

open('right_ligands.txt', 'w').writelines(smiles + '\n' for smiles in right_ligands)
