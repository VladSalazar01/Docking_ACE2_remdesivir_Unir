#!/usr/bin/env python3
"""
Calcula el centro del grid desde crystal_ligand.mol2
Carpeta: G:\Cribado_VC_ACE_2\validacion_ACE1
"""

from pathlib import Path

WORK_DIR = Path(r"G:\Cribado_VC_ACE_2\validacion_ACE1")
LIGAND_FILE = WORK_DIR / "crystal_ligand.mol2"

def calcular_centroide(mol2_file):
    coords = []
    in_atom_section = False
    with open(mol2_file, 'r') as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atom_section = True
                continue
            if '@<TRIPOS>' in line and in_atom_section:
                break
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                        coords.append((x, y, z))
                    except ValueError:
                        continue
    if not coords:
        print("ERROR: No se encontraron coordenadas")
        return None
    cx = sum(c[0] for c in coords) / len(coords)
    cy = sum(c[1] for c in coords) / len(coords)
    cz = sum(c[2] for c in coords) / len(coords)
    return cx, cy, cz, len(coords)

if __name__ == "__main__":
    print("="*50)
    print("CALCULO DE CENTRO DEL GRID")
    print("="*50)
    result = calcular_centroide(LIGAND_FILE)
    if result:
        cx, cy, cz, n_atoms = result
        print(f"\nLigando: {LIGAND_FILE.name}")
        print(f"Atomos: {n_atoms}")
        print(f"\nCENTRO DEL GRID:")
        print(f"  X = {cx:.3f}")
        print(f"  Y = {cy:.3f}")
        print(f"  Z = {cz:.3f}")
        print(f"\nActualizar en run_gnina_validation.py:")
        print(f"  CENTER = ({cx:.1f}, {cy:.1f}, {cz:.1f})")