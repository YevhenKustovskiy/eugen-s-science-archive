import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser(description="Assigns bond records based on template structure using RDKit")
parser.add_argument("-t", type=str, help="Specify template file")
parser.add_argument("-f", type=str, help="Specify file, which contains conformation of molecule to assign bonds to")
parser.add_argument("-o", type=str, help="Specify output filename")
args = parser.parse_args()

if __name__ == "__main__":
    template_mol = Chem.MolFromMolFile(args.t)
    
    if template_mol is None:
        raise ValueError("Failed to load template molecule. Check file format and path.")
        
    molecule = Chem.MolFromMolFile(args.f)
    
    if molecule is None:
        raise ValueError("Failed to load molecule. Check file format and path.")
    
    new_molecule = AllChem.AssignBondOrdersFromTemplate(template_mol, molecule)
    Chem.MolToMolFile(new_molecule, args.o)
    
