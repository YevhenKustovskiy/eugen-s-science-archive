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
                
        