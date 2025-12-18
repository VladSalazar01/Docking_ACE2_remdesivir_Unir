
from pathlib import Path
import numpy as np

cif = Path(r'G:\Boltz_benchmark\boltz_outputs\run_001\boltz_results_input_mini\predictions\input_mini\input_mini_model_0.cif')
ref = Path(r'G:\Boltz_benchmark\ligand_A1IDX.pdb')

# Extraer ligando predicho
coords_pred = []
for line in cif.read_text().splitlines():
    if line.startswith('HETATM') and 'LIG1' in line:
        parts = line.split()
        if len(parts) >= 13 and not parts[3].startswith('H'):
            coords_pred.append([float(parts[10]), float(parts[11]), float(parts[12])])

# Extraer ligando referencia
coords_ref = []
for line in ref.read_text().splitlines():
    if line.startswith(('ATOM', 'HETATM')) and not line[12:16].strip().startswith('H'):
        coords_ref.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])

print(f'Referencia: {len(coords_ref)} átomos')
print(f'Predicho: {len(coords_pred)} átomos')
