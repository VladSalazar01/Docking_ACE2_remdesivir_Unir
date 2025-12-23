#!/usr/bin/env python3
"""Generate DOCK6 spheres centered on binding site - No numpy needed"""

# Centro del sitio de unión (de Gnina Fragment 1)
center_x, center_y, center_z = 16.48, 15.24, 25.54
box_size = 20.0  # Angstroms
spacing = 2.0    # Espaciado entre esferas
radius = 1.8     # Radio de cada esfera

spheres = []
half_box = box_size / 2

# Generar posiciones
x = -half_box
while x <= half_box:
    y = -half_box
    while y <= half_box:
        z = -half_box
        while z <= half_box:
            px = center_x + x
            py = center_y + y
            pz = center_z + z
            spheres.append((px, py, pz))
            z += spacing
        y += spacing
    x += spacing

# Escribir archivo de esferas
with open('selected_spheres.sph', 'w') as f:
    f.write('DOCK spheres centered on Fragment 1 binding site\n')
    f.write(f'cluster     1   number of spheres in cluster   {len(spheres)}\n')
    for i, (sx, sy, sz) in enumerate(spheres, 1):
        f.write(f'{i:6d}  {sx:8.3f}  {sy:8.3f}  {sz:8.3f}  {radius:.3f}  0  0  0\n')

print(f"Generated {len(spheres)} spheres")
print(f"Center: ({center_x}, {center_y}, {center_z})")
print(f"Box size: {box_size} A")
