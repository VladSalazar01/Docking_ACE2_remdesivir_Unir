#!/usr/bin/env python3
"""
Parsear CIF en vez de PDB (9FMM no tiene PDB disponible)
"""
from pathlib import Path

def parse_cif_ligands():
    base = Path(r"G:\Smina_benchmark")
    cif_file = base / "9FMM_original.cif"
    
    if not cif_file.exists():
        print("CIF not found")
        return
    
    print(f"Parsing {cif_file.name}...\n")
    
    # Encontrar ligandos no estándar
    ligands = {}
    in_atom_site = False
    
    with open(cif_file) as f:
        for line in f:
            if line.startswith('_atom_site.'):
                in_atom_site = True
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) > 5:
                    atom_type = parts[0]
                    comp_id = parts[5]  # Residue name
                    
                    if atom_type == 'HETATM' and comp_id not in ['HOH', 'WAT']:
                        if comp_id not in ligands:
                            ligands[comp_id] = 0
                        # Contar solo heavy atoms
                        element = parts[-1] if len(parts) > 10 else parts[2]
                        if element != 'H':
                            ligands[comp_id] += 1
    
    print(f"{'='*60}")
    print("LIGANDS IN CRYSTAL:")
    print(f"{'='*60}")
    
    for lig, count in ligands.items():
        print(f"{lig}: {count} heavy atoms")
        
        # Comparar con tu fragment2 (28 heavy atoms)
        if count == 28:
            print(f"  ✓ MATCHES your fragment2!")
        elif count > 28:
            print(f"  (Your fragment2 has {28} atoms - this is larger)")
        else:
            print(f"  (Your fragment2 has {28} atoms - this is smaller)")

if __name__ == "__main__":
    parse_cif_ligands()