#!/usr/bin/env python3
"""Generate DOCK6 spheres from MSMS surface, filtered by binding site"""

CENTER = (16.48, 15.24, 25.54)
SITE_RADIUS = 12.0
SPHERE_RADIUS = 1.4
MIN_SPACING = 1.5

def distance(p1, p2):
    return ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)**0.5

vertices = []
with open('receptor_surface.vert', 'r') as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 3:
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                vertices.append((x, y, z))
            except ValueError:
                continue

print(f"Total surface vertices: {len(vertices)}")

site_vertices = [v for v in vertices if distance(v, CENTER) <= SITE_RADIUS]
print(f"Vertices near binding site: {len(site_vertices)}")

spheres = []
for v in site_vertices:
    too_close = False
    for s in spheres:
        if distance(v, s) < MIN_SPACING:
            too_close = True
            break
    if not too_close:
        spheres.append(v)

print(f"Final spheres: {len(spheres)}")

with open('selected_spheres.sph', 'w') as f:
    f.write('DOCK spheres from MSMS surface\n')
    f.write(f'cluster     1   number of spheres in cluster   {len(spheres)}\n')
    for i, (x, y, z) in enumerate(spheres, 1):
        f.write(f'{i:6d}  {x:8.3f}  {y:8.3f}  {z:8.3f}  {SPHERE_RADIUS:.3f}  0  0  0\n')

print(f"Saved to selected_spheres.sph")
