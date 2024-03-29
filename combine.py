from rdkit import Chem
from fragment import fragment
from ligate import ligate

left_ligands = [Chem.MolFromSmiles(line.strip()) for line in open('left_ligands.txt')]
right_ligands = [Chem.MolFromSmiles(line.strip()) for line in open('right_ligands.txt')]

with open('tridentate_carbene_smiles.txt', 'w') as outfile:
    for left_ligand in left_ligands:
        for right_ligand in right_ligands:
            left_ligand_noir = fragment(left_ligand)[0]
            right_ligand_noir = fragment(right_ligand)[0]
            molecule = ligate([left_ligand_noir, right_ligand_noir])
            molecule_smiles = Chem.MolToSmiles(molecule)
            outfile.write(molecule_smiles + '\n')
