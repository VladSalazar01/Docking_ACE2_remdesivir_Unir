#!/usr/bin/env python3
"""Calcula RMSD entre pose predicha y referencia"""

def read_mol2_coords(filename):
    coords = []
    in_atoms = False
    with open(filename, 'r') as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atoms = True
                continue
            elif '@<TRIPOS>' in line:
                in_atoms = False
                continue
            if in_atoms:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        atom_type = parts[5] if len(parts) > 5 else parts[1]
                        # Solo atomos pesados (no H)
                        if not atom_type.startswith('H'):
                            coords.append((x, y, z))
                    except:
                        pass
    return coords

def read_pdb_coords(filename):
    coords = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                if not atom_name.startswith('H'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
    return coords

# Leer coordenadas
ref_coords = read_pdb_coords('ligand_ref.pdb')
pred_coords = read_mol2_coords('output_run1_scored.mol2')

print(f"Atomos pesados referencia: {len(ref_coords)}")
print(f"Atomos pesados prediccion: {len(pred_coords)}")

if len(ref_coords) == len(pred_coords):
    # Calcular RMSD
    sum_sq = 0
    for (rx, ry, rz), (px, py, pz) in zip(ref_coords, pred_coords):
        sum_sq += (rx-px)**2 + (ry-py)**2 + (rz-pz)**2
    rmsd = (sum_sq / len(ref_coords)) ** 0.5
    print(f"\nRMSD: {rmsd:.3f} A")
else:
    print("ERROR: Numero de atomos diferente!")
