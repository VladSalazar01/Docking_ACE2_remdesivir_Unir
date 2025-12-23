#!/usr/bin/env python3
"""Extrae residuos del receptor dentro de un radio del ligando"""

# Leer coordenadas del ligando
lig_coords = []
with open('ligand_ref.pdb', 'r') as f:
    for line in f:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            lig_coords.append((x, y, z))

print(f"Ligando: {len(lig_coords)} atomos")

# Calcular centro del ligando
lig_cx = sum(c[0] for c in lig_coords) / len(lig_coords)
lig_cy = sum(c[1] for c in lig_coords) / len(lig_coords)
lig_cz = sum(c[2] for c in lig_coords) / len(lig_coords)
print(f"Centro ligando: ({lig_cx:.2f}, {lig_cy:.2f}, {lig_cz:.2f})")

CUTOFF = 15.0  # Angstroms desde cualquier atomo del ligando

def min_dist_to_ligand(x, y, z):
    """Distancia minima a cualquier atomo del ligando"""
    min_d = 9999
    for lx, ly, lz in lig_coords:
        d = ((x-lx)**2 + (y-ly)**2 + (z-lz)**2)**0.5
        if d < min_d:
            min_d = d
    return min_d

# Leer receptor mol2 y filtrar
in_atoms = False
in_bonds = False
atoms_kept = []
atom_indices_kept = set()
all_lines = []
header_lines = []
atom_lines = []
bond_lines = []
other_lines = []

with open('receptor.mol2', 'r') as f:
    current_section = 'header'
    for line in f:
        if '@<TRIPOS>ATOM' in line:
            current_section = 'atoms'
            continue
        elif '@<TRIPOS>BOND' in line:
            current_section = 'bonds'
            continue
        elif line.startswith('@<TRIPOS>'):
            current_section = 'other'
            other_lines.append(line)
            continue
        
        if current_section == 'header':
            header_lines.append(line)
        elif current_section == 'atoms':
            atom_lines.append(line)
        elif current_section == 'bonds':
            bond_lines.append(line)
        else:
            other_lines.append(line)

print(f"Receptor original: {len(atom_lines)} atomos")

# Filtrar atomos cercanos al ligando
kept_atoms = []
old_to_new_idx = {}
new_idx = 1

for line in atom_lines:
    parts = line.split()
    if len(parts) >= 6:
        try:
            old_idx = int(parts[0])
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            
            if min_dist_to_ligand(x, y, z) <= CUTOFF:
                old_to_new_idx[old_idx] = new_idx
                # Reescribir con nuevo indice
                parts[0] = str(new_idx)
                new_line = f"{new_idx:7d} {parts[1]:<4s}    {x:9.4f} {y:9.4f} {z:9.4f} {' '.join(parts[5:])}\n"
                kept_atoms.append(new_line)
                new_idx += 1
        except:
            pass

print(f"Receptor recortado: {len(kept_atoms)} atomos (dentro de {CUTOFF} A del ligando)")

# Filtrar bonds que conectan atomos mantenidos
kept_bonds = []
bond_idx = 1
for line in bond_lines:
    parts = line.split()
    if len(parts) >= 4:
        try:
            a1 = int(parts[1])
            a2 = int(parts[2])
            if a1 in old_to_new_idx and a2 in old_to_new_idx:
                new_a1 = old_to_new_idx[a1]
                new_a2 = old_to_new_idx[a2]
                new_line = f"{bond_idx:6d} {new_a1:5d} {new_a2:5d} {parts[3]}\n"
                kept_bonds.append(new_line)
                bond_idx += 1
        except:
            pass

print(f"Bonds mantenidos: {len(kept_bonds)}")

# Escribir nuevo archivo
with open('receptor_site.mol2', 'w') as f:
    f.write('@<TRIPOS>MOLECULE\n')
    f.write('receptor_binding_site\n')
    f.write(f'{len(kept_atoms)} {len(kept_bonds)} 0 0 0\n')
    f.write('SMALL\n')
    f.write('GASTEIGER\n')
    f.write('\n')
    f.write('@<TRIPOS>ATOM\n')
    for line in kept_atoms:
        f.write(line)
    f.write('@<TRIPOS>BOND\n')
    for line in kept_bonds:
        f.write(line)

print("Guardado: receptor_site.mol2")
