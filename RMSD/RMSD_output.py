#Designed to be run in PyMol with the following structures loaded
from pymol import cmd
from itertools import combinations

# Load the structures
structures = ['e_faecalis', 'human_zinc', 'MRSA', 's_pneumoniae', 's_aureus']
for structure in structures:
    cmd.load(f'{structure}.pdb')
    print(f"Loaded {structure}.")

# Open a file to write the RMSD values
with open('rmsd_output.txt', 'w') as f:
    # Calculate and write RMSD for each pair of structures
    for structure1, structure2 in combinations(structures, 2):
        rmsd = cmd.align(f'{structure1}', f'{structure2}')[0]
        f.write(f"RMSD between {structure1} and {structure2}: {rmsd}\n")
        print(f"Calculated RMSD between {structure1} and {structure2}.")

