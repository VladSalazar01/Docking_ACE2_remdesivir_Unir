#!/usr/bin/env python3
"""
Debug: Comparar coordenadas extraídas
"""
from pathlib import Path
import numpy as np

def extract_coords(pdb_file):
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ENDMDL'):
                break
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                if atom_name.startswith('H'):
                    continue
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except:
                    continue
    return np.array(coords)

def simple_rmsd(c1, c2):
    if len(c1) != len(c2):
        return 999.0
    diff = c1 - c2
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))

base = Path(r"G:\Smina_benchmark")

# Referencia
ref = base / "ligand_A1IDX_fragment2.pdb"
ref_coords = extract_coords(ref)

print(f"Reference: {ref}")
print(f"Atoms: {len(ref_coords)}")
print(f"Center: {ref_coords.mean(axis=0)}")
print(f"Coords range:")
print(f"  X: {ref_coords[:,0].min():.2f} to {ref_coords[:,0].max():.2f}")
print(f"  Y: {ref_coords[:,1].min():.2f} to {ref_coords[:,1].max():.2f}")
print(f"  Z: {ref_coords[:,2].min():.2f} to {ref_coords[:,2].max():.2f}")

# Mejor docking
docked = base / "grid_results_seq/config_003/docked.pdb"
if docked.exists():
    docked_coords = extract_coords(docked)
    
    print(f"\nDocked (best): {docked}")
    print(f"Atoms: {len(docked_coords)}")
    print(f"Center: {docked_coords.mean(axis=0)}")
    print(f"Coords range:")
    print(f"  X: {docked_coords[:,0].min():.2f} to {docked_coords[:,0].max():.2f}")
    print(f"  Y: {docked_coords[:,1].min():.2f} to {docked_coords[:,1].max():.2f}")
    print(f"  Z: {docked_coords[:,2].min():.2f} to {docked_coords[:,2].max():.2f}")
    
    # RMSD sin alineación
    rmsd_raw = simple_rmsd(ref_coords, docked_coords)
    print(f"\nRMSD (no alignment): {rmsd_raw:.3f} Å")
    
    # Centrar ambos en origen
    ref_centered = ref_coords - ref_coords.mean(axis=0)
    docked_centered = docked_coords - docked_coords.mean(axis=0)
    
    rmsd_centered = simple_rmsd(ref_centered, docked_centered)
    print(f"RMSD (centered): {rmsd_centered:.3f} Å")