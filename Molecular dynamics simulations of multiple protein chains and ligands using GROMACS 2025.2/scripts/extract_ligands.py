"""extract_ligands.py: separates ligand coordinates from complex coordinates
   Copyright (C) 2025  Yevhen Kustovskiy
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 
   
   Copies of GNU license, NumPy, RDKit, and MDAnalysis license are available at:
   
   https://github.com/YevhenKustovskiy/eugen-s-science-archive/blob/main/LICENSE.txt 
   
   Contact email: ykustovskiy@gmail.com
   
   Requirements: script was succesfully used with 
   Python 3.11, NumPy 1.24.3, RDKit 2023.9.4, MDAnalysis 2.5.0, 
"""

import os
import shutil
import argparse
import warnings
import numpy as np
import MDAnalysis as mda

from rdkit.Chem import MolToMolFile


parser = argparse.ArgumentParser(description = "Extracts system(s) components into separate files")
parser.add_argument("-p",  help = "Define complex file name ", type=str)
parser.add_argument("-batch",  help = "Extract components of all files in directory",
                    action="store_true")
parser.add_argument("-nowarn", help="Mute warnings which may occure during evaluation", action="store_true")
parser.add_argument("-v", help="Provide the information about progress", action="store_true")
args = parser.parse_args()

fnames = []
if args.batch:
    for file in os.listdir():
        if "pdb" in file:
            fnames.append(file)

elif args.p:
    fnames.append(args.p)

def extract(fname: str):
    system = mda.Universe(fname)

    outdir = f"{fname.split('.')[0]}_components"
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.mkdir(outdir)

    protein = system.select_atoms("protein")
    protein.write(os.path.join(outdir,"protein.pdb"))

    hetero = system.atoms - protein
    hetero_names = np.unique(hetero.residues.resnames)

    for name in hetero_names:
        res = hetero.select_atoms(f"resname {name}")
        mol = res.atoms.convert_to("RDKIT", force=True)
        MolToMolFile(mol, os.path.join(outdir,f"{name}.mol"))


if __name__ == "__main__":
    if args.nowarn:
        warnings.filterwarnings("ignore")
        print("Warnings are muted!")

    for file in fnames:
        if args.v:
            print(f"Processing file {file} ...")
        extract(file)










