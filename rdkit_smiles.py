#!/usr/bin/env python
"""
Script to read in Smiles strings with RDkit.
Author: Richard Bradshaw, richard.bradshaw@nih.gov
"""

# AllChem has 2D/3D conformer generation
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from pandas import read_csv

mol_file = "/u/bradshaw/WORK/Git/SAMPL6/physical_properties/pKa/molecule_ID_and_SMILES.csv"

mol_table = read_csv(mol_file)
mols = []

for idx, name, smiles in mol_table.itertuples():
    tmpmol = Chem.MolFromSmiles(smiles)
    tmpmol.SetProp("_Name", name)
    #tmpmol.GetPropsAsDict(includePrivate=True) # Dict of defined molecule properties (currently only _Name!)
    tmpmol = Chem.AddHs(tmpmol)
    print(smiles)
    print('Total atoms for molecule %s = %d' % (name, tmpmol.GetNumAtoms()))
    mols.append(tmpmol)

# Create 2D depictions & 3D models
for currmol in mols:
    Chem.Compute2DCoords(currmol)
    Draw.MolToFile(currmol, "RDkit_%s_2D.png" % currmol.GetProp("_Name"))
    Chem.EmbedMolecule(currmol)
    Chem.MMFFOptimizeMolecule(currmol)
    Chem.MolToPDBFile(currmol, "RDkit_%s_3D.pdb" % currmol.GetProp("_Name"))


