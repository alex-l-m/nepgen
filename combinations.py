import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Load substitutions from a file
substitutions = pd.read_csv('substitutions.csv')
# Column 'smiles' as a list
substitution_smiles = substitutions['smiles'].tolist()
sub_1 = substitution_smiles
sub_3 = substitution_smiles

ligand_l = Chem.MolFromSmiles('SC1=C(S)N=C(N2C([Ir]3C4N(C5=NC(S)=C(S)N=C5N4C6=C(C(I)=C(C2=C63)Cl)Cl)C)N7C)C7=N1')


def sub_combinations(ligand, sub1, sub2, sub3,              # returns a set of all combinations of substitutions
                     p1=Chem.MolFromSmiles('S'),            # limitation, 3 lists, 3 groups to substitute
                     p2=Chem.MolFromSmiles('Cl'),
                     p3=Chem.MolFromSmiles('I')):
    substituted_ligands = set()
    for x in sub1:
        for y in sub2:
            for z in sub3:
                ligand_w_sub1 = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
                ligand_w_sub2 = AllChem.ReplaceSubstructs(ligand_w_sub1, p2, Chem.MolFromSmiles(y), True)[0]
                ligand_w_sub3 = AllChem.ReplaceSubstructs(ligand_w_sub2, p3, Chem.MolFromSmiles(z), True)[0]
                if ligand_w_sub3 is not None:
                    substituted_ligands.add(Chem.MolToSmiles(ligand_w_sub3, True))
    return substituted_ligands


sub_2 = ["[CH3]", "[H]"]

set1_unfiltered = sub_combinations(ligand_l, sub_1, sub_2, sub_3)
set1 = set(smiles for smiles in set1_unfiltered if Chem.MolFromSmiles(smiles) is not None)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set1], molsPerRow=4, subImgSize=(200, 200)).show()

# Save SMILES strings of substituted ligands
output_smiles_filename = 'substituted_ligands.txt'
with open(output_smiles_filename, 'w') as f:
    for x in set1:
        f.write(x + '\n')
