"""add_structure.py: merges two coordinate files
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
   
   Copies of GNU license and MDAnalysis license are available at:
   
   https://github.com/YevhenKustovskiy/eugen-s-science-archive/blob/main/LICENSE.txt 
   
   Contact email: ykustovskiy@gmail.com
   
   Requirements: script was succesfully used with 
   Python 3.10, MDAnalysis 2.5.0
"""

import argparse
import MDAnalysis as mda


parser = argparse.ArgumentParser(description = 'Appends coordinates from one file to another!')
parser.add_argument("-p",  help = "Define protein structure file name ")
parser.add_argument("-l",  help = "Define ligand structure file name ")
parser.add_argument("-o",  help = "Specify the output file name.", default = "out.pdb")
args = parser.parse_args()


prot = mda.Universe(args.p)
lig  = mda.Universe(args.l)
sys  = mda.Merge(prot.atoms, lig.atoms)
sys.atoms.write(args.o)
