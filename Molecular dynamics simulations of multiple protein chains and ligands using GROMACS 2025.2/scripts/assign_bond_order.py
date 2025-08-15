"""assign_bond_order: translates bonds from template molecule to different conformation of the same molecule
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
   
   A copies of GNU license, RDKit license are available at:
   
   https://github.com/YevhenKustovskiy/md-scripts/edit/main/LICENSE   
   
   Contact email: ykustovskiy@gmail.com
   
   Requirements: script was succesfully implemented with 
   Python 3.10, RDKit 2023.9.4 
"""

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
    

