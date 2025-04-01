import sys
import os
import re

# Mapping of three-letter amino acid codes to one-letter codes
aa3to1 = {
    'ALA': 'A', 'VAL': 'V', 'PHE': 'F', 'PRO': 'P', 'MET': 'M',
    'ILE': 'I', 'LEU': 'L', 'ASP': 'D', 'GLU': 'E', 'LYS': 'K',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'HIS': 'H',
    'CYS': 'C', 'ASN': 'N', 'GLN': 'Q', 'TRP': 'W', 'GLY': 'G',
    'MSE': 'M',  # Selenomethionine treated as Methionine
}

# Regular expression pattern to match CA atoms in PDB files
ca_pattern = re.compile(r"^(ATOM|HETATM)\s+\d+\s+CA\s+[A-Z]?\s*([A-Z]{3})\s([\s\w])")

def parse_pdb(pdb_path):
    """
    Parses a PDB file and extracts chain sequences.
    """
    chain_dict = {}
    chain_list = []
    try:
        with open(pdb_path, 'r') as fp:
            for line in fp:
                if line.startswith("ENDMDL"):  # Stop processing after the first model
                    break
                match = ca_pattern.match(line)
                if match:
                    resn, chain = match.group(2), match.group(3).strip()
                    if resn in aa3to1:
                        if chain in chain_dict:
                            chain_dict[chain] += aa3to1[resn]
                        else:
                            chain_dict[chain] = aa3to1[resn]
                            chain_list.append(chain)
        return chain_dict, os.path.basename(pdb_path).split('.')[0]
    except FileNotFoundError:
        sys.stderr.write(f"Error: File {pdb_path} not found.\n")
        return None, None
    except Exception as e:
        sys.stderr.write(f"Error processing {pdb_path}: {e}\n")
        return None, None

def pdb_to_fasta(pdb_file, output_file=None):
    """Converts a single PDB file into FASTA format."""
    chain_dict, fasta_header = parse_pdb(pdb_file)
    if chain_dict is None:
        return
    
    if output_file:
        with open(output_file, 'a') as output:
            for chain, sequence in chain_dict.items():
                output.write(f">{fasta_header}:{chain}\n{sequence}\n")
    else:
        for chain, sequence in chain_dict.items():
            sys.stdout.write(f">{fasta_header}:{chain}\n{sequence}\n")

def pdbs_to_fasta(folder_path, output_file):
    """Processes all PDB files in a folder and writes the FASTA sequences to an output file."""
    try:
        with open(output_file, 'w') as output:
            for pdb_file in os.listdir(folder_path):
                if pdb_file.endswith(".pdb"):
                    pdb_to_fasta(os.path.join(folder_path, pdb_file), output_file)
        print(f"FASTA sequences saved to {output_file}")
    except Exception as e:
        sys.stderr.write(f"Error processing folder {folder_path}: {e}\n")

if __name__ == "__main__":
    if sys.argv[1] == '-h': # Help message
        sys.stderr.write("Usage:\n")
        sys.stderr.write("  Convert a single PDB file to FASTA: pdb2fasta.py [input.pdb] [output.fasta]\n")
        sys.stderr.write("  Convert all PDB files in a folder: pdb2fasta.py --folder [input_folder] [output.fasta]\n")
        sys.exit(1)
    
    elif sys.argv[1] == "--folder" and len(sys.argv) == 4: # Process all PDB files in a folder
        sys.stdout.write(f"Processing all PDB file in folder {sys.argv[2]} \n")
        pdbs_to_fasta(sys.argv[2], sys.argv[3])
    else: # Process a single PDB file
        sys.stdout.write(f"Processing single PDB file {sys.argv[1]} \n")
        pdb_to_fasta(sys.argv[1], sys.argv[2])
