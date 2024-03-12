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


set1_unfiltered = sub_combinations(ligand_l, sub_1, sub_2, sub_3)
set1 = set(smiles for smiles in set1_unfiltered if Chem.MolFromSmiles(smiles) is not None)
Draw.MolsToGridImage([Chem.AddHs(Chem.MolFromSmiles(x)) for x in set1], molsPerRow=8, subImgSize=(200, 200)).show()


def top_combinations(ligand, top_sub, p1=Chem.MolFromSmiles('S')):
    top_ligands = set()
    for x in top_sub:
        ligand_w_top_sub = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
        if ligand_w_top_sub is not None:
            top_ligands.add(Chem.MolToSmiles(ligand_w_top_sub, True))
    return top_ligands


def bot_combinations(ligand, bot_sub, p1=Chem.MolFromSmiles('S')):
    bot_ligands = set()
    for x in bot_sub:
        ligand_w_bot_sub = AllChem.ReplaceSubstructs(ligand, p1, Chem.MolFromSmiles(x), True)[0]
        if ligand_w_bot_sub is not None:
            bot_ligands.add(Chem.MolToSmiles(ligand_w_bot_sub, True))
    return bot_ligands


# # Save SMILES strings of substituted ligands
# output_smiles_filename = 'substituted_ligands.txt'
# with open(output_smiles_filename, 'w') as f:
#     for x in set1:
#         f.write(x + '\n')
