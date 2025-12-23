#!/usr/bin/env python3
"""Convert PDB to XYZR with standard atomic radii"""

# Radios atómicos estándar (van der Waals)
RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
    'H': 1.20, 'P': 1.80, 'F': 1.47, 'CL': 1.75,
    'BR': 1.85, 'I': 1.98, 'ZN': 1.39, 'FE': 1.40,
    'CA': 1.97, 'MG': 1.73, 'NA': 2.27, 'K': 2.75,
}

with open('receptor.pdb', 'r') as f_in, open('receptor.xyzr', 'w') as f_out:
    for line in f_in:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            
            # Obtener elemento del PDB
            element = line[76:78].strip().upper()
            if not element:
                # Inferir del nombre del átomo
                atom_name = line[12:16].strip()
                element = atom_name[0] if atom_name[0] in 'CNOSPFH' else 'C'
            
            radius = RADII.get(element, 1.70)
            f_out.write(f"{x:8.3f} {y:8.3f} {z:8.3f} {radius:.2f}\n")

print("XYZR file generated with correct radii")
