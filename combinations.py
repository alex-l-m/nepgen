import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

# Load substitutions from a file
substitutions = pd.read_csv('substitutions.csv')
# Column 'smiles' as a list
substitution_smiles = substitutions['smiles'].tolist()
sub_1 = substitution_smiles
sub_2 = ["[CH3]", "[H]"]
sub_3 = substitution_smiles

ligand_l = Chem.MolFromSmiles('SC1=C(S)N=C(N2C([Ir]3C4N(C5=NC(S)=C(S)N=C5N4C6=C(C(I)=C(C2=C63)Cl)Cl)C)N7C)C7=N1')
ligand_top = Chem.MolFromSmiles('SC1=CN=C(C2=N1)N3C(N2C)[Ir]4C5N(C)C6=NC(S)=CN=C6N5C7=CC(C(F)(F)F)=CC3=C74')
ligand_bot = Chem.MolFromSmiles('SC1=CN=C2C(N3C(N2C)[Ir]4C5N(C)C6=NC=C(N=C6N5C7=CC(C(F)(F)F)=CC3=C74)S)=N1')


def smiles_to_rdkit_mol(sub_list):  # converts list of smiles strings to rdkit molecules
    rdkit_list = []                 # list of rdkit molecules
    for sub in sub_list:            # goes through lists and does conversion
        rdkit_mol = Chem.MolFromSmiles(sub)
        rdkit_list.append(rdkit_mol)
    return rdkit_list


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


def sub_one_group(ligand, sub, p1=Chem.MolFromSmiles('S')):
    substituted_ligands = set()
    for x in sub:
        ligand_w_sub = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
        if ligand_w_sub is not None:
            substituted_ligands.add(Chem.MolToSmiles(ligand_w_sub, True))
    return substituted_ligands


# # Save SMILES strings of substituted ligands
# output_smiles_filename = 'substituted_ligands.txt'
# with open(output_smiles_filename, 'w') as f:
#     for x in set1:
#         f.write(x + '\n')

# Show images of sub_combinations
# set1_unfiltered = sub_combinations(ligand_l, sub_1, sub_2, sub_3)
# set1 = set(smiles for smiles in set1_unfiltered if Chem.MolFromSmiles(smiles) is not None)
# Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set1], molsPerRow=8, subImgSize=(200, 200)).show()

# Show images of sub_one_group, issues with mesityl group
sub_one_set = sub_one_group(ligand_l, sub_1)
sos_set = set(smiles for smiles in sub_one_set if Chem.MolFromSmiles(smiles) is not None)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in sos_set], molsPerRow=8, subImgSize=(200, 200)).show()
