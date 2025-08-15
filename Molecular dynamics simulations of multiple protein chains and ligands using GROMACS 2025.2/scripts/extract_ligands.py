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







