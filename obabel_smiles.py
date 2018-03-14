#!/usr/bin/env python
"""
Script to read in Smiles strings with openbabel.
Author: Richard Bradshaw, richard.bradshaw@nih.gov
"""

import openbabel as ob
import pybel as pyb
from pandas import read_csv

mol_file = "/u/bradshaw/WORK/Git/SAMPL6/physical_properties/pKa/molecule_ID_and_SMILES.csv"

mol_table = read_csv(mol_file)
mols = []
obConverter = ob.OBConversion()
obConverter.SetInAndOutFormats("smi", "pdb")

for idx, name, smiles in mol_table.itertuples():
    tmpmol = ob.OBMol()
    obConverter.ReadString(tmpmol, smiles)
    tmpmol.AddHydrogens()
    print(smiles)
    print('Total atoms for molecule %s = %d' % (name, tmpmol.NumAtoms()))
    mols.append(tmpmol)
    obConverter.WriteFile(tmpmol, '%s.pdb' % name)

pybmols = map(pyb.Molecule,mols)
map(lambda x: x.make3D(),pybmols)
for curr, name in zip(pybmols, mol_table['SAMPL6 Molecule ID']):
    curr.title = name
    curr.write('pdb', '%s_3D.pdb' % name, overwrite=True)
    curr.draw(show=False, filename='%s_2D.png' % name)


