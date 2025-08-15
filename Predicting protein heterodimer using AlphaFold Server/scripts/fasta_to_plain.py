"""fasta_to_plain: in batch converts fasta files with aminoacid sequences to lines with no symbols other than aminoacids 
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
   
   A copy of GNU license is available at:
   
   https://github.com/YevhenKustovskiy/eugen-s-science-archive/edit/main/LICENSE   
   
   Contact email: ykustovskiy@gmail.com
   
   Requirements: Python 3.11
"""

import os
import argparse

# example of input (Windows): fast_to_plain.py -f r"C:\Absolute\Path\To\File.fasta"
# example of input (Windows): fast_to_plain.py -d
# example of input (Linux/Ubuntu): python3 fast_to_plain.py -f /absolute/path/to/file.fasta
# example of input (Linux/Ubuntu): python3 fast_to_plain.py -d

CWD = os.getcwd();
EXT = ".fasta";

def process_file(filename):
    with open(filename, 'r', encoding='utf-8') as infile, \
         open(os.path.splitext(filename)[0] + ".plain", 'w', encoding='utf-8') as outfile:
        
        next(infile)
        outfile.write(''.join(line.strip() for line in infile))

def process_folder(folder):
    for file in os.listdir(CWD):
        if os.path.isfile(file):
            if not file.endswith(EXT):
                continue
            
            process_file(file)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts .fasta file into plain text file with no technical symbols (e.g.)")
    
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("-f", help="Define file (.fasta)")
    group.add_argument("-d", help="Process all files in directory", action="store_true")
    
    args = parser.parse_args()
    
    if args.f:
        process_file(args.f)
    elif args.d:
        process_folder(CWD)
                

        


