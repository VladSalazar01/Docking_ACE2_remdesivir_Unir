#!/usr/bin/env python3
"""Genera esferas DOCK6 basadas en las posiciones de los átomos del ligando"""

# Leer coordenadas del ligando PDB
coords = []
with open('ligand_ref.pdb', 'r') as f:
    for line in f:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append((x, y, z))

print(f"Leidos {len(coords)} atomos del ligando")

# Calcular centro del ligando
cx = sum(c[0] for c in coords) / len(coords)
cy = sum(c[1] for c in coords) / len(coords)
cz = sum(c[2] for c in coords) / len(coords)
print(f"Centro del ligando: ({cx:.3f}, {cy:.3f}, {cz:.3f})")

# Generar esferas: una por cada atomo + algunas adicionales en el centro
spheres = []
radius = 1.8

# Esferas en posiciones atomicas
for x, y, z in coords:
    spheres.append((x, y, z, radius))

# Agregar esferas adicionales alrededor del centro (para mejor cobertura)
offsets = [0, -2, 2, -4, 4]
for dx in offsets:
    for dy in offsets:
        for dz in offsets:
            if abs(dx) + abs(dy) + abs(dz) <= 6:  # Limitar a esfera aproximada
                spheres.append((cx + dx, cy + dy, cz + dz, radius))

# Eliminar duplicados cercanos
unique_spheres = []
for s in spheres:
    is_dup = False
    for us in unique_spheres:
        dist = ((s[0]-us[0])**2 + (s[1]-us[1])**2 + (s[2]-us[2])**2)**0.5
        if dist < 1.5:
            is_dup = True
            break
    if not is_dup:
        unique_spheres.append(s)

print(f"Generadas {len(unique_spheres)} esferas unicas")

# Escribir archivo
with open('selected_spheres.sph', 'w') as f:
    f.write('DOCK spheres based on ligand atom positions\n')
    f.write(f'cluster     1   number of spheres in cluster   {len(unique_spheres)}\n')
    for i, (x, y, z, r) in enumerate(unique_spheres, 1):
        f.write(f'{i:6d}  {x:8.3f}  {y:8.3f}  {z:8.3f}  {r:.3f}  0  0  0\n')

print(f"Archivo guardado: selected_spheres.sph")
