#!/usr/bin/env python3
"""
Cálculo correcto de RMSD para ligando completo
"""
from pathlib import Path
import numpy as np

def extract_coords(pdb):
    coords = []
    with open(pdb, 'r') as f:
        for line in f:
            if line.startswith('ENDMDL'):
                break
            if line.startswith(('ATOM', 'HETATM')) and not line[12:16].strip().startswith('H'):
                try:
                    coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except:
                    pass
    return np.array(coords)

def kabsch_rmsd(P, Q):
    if len(P) != len(Q):
        return 999.0
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    C = np.dot(Q_c.T, P_c)
    V, S, W = np.linalg.svd(C)
    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]
    Q_rot = np.dot(Q_c, np.dot(V, W))
    return np.sqrt(np.mean(np.sum((P_c - Q_rot)**2, axis=1)))

base = Path(r"G:\Smina_benchmark")

# Fragment2
frag2_ref = extract_coords(base / "ligand_A1IDX_fragment2.pdb")
frag2_docked = extract_coords(base / "multiseed/run_000/docked.pdb")
rmsd_frag2 = kabsch_rmsd(frag2_ref, frag2_docked)

# Complete
complete_ref = extract_coords(base / "ligand_A1IDX.pdb")
complete_docked = extract_coords(base / "debug_complete/docked.pdb")
rmsd_complete = kabsch_rmsd(complete_ref, complete_docked)

print("="*60)
print("RMSD COMPARISON: Fragment2 vs Complete")
print("="*60)
print(f"\nFragment2:")
print(f"  Reference atoms: {len(frag2_ref)}")
print(f"  Docked atoms: {len(frag2_docked)}")
print(f"  RMSD: {rmsd_frag2:.3f} Å")

print(f"\nComplete:")
print(f"  Reference atoms: {len(complete_ref)}")
print(f"  Docked atoms: {len(complete_docked)}")
print(f"  RMSD: {rmsd_complete:.3f} Å")

print(f"\n{'='*60}")
if rmsd_complete < rmsd_frag2:
    print(f"✓ Complete ligand is BETTER by {rmsd_frag2 - rmsd_complete:.3f} Å")
else:
    print(f"✗ Fragment2 remains better by {rmsd_complete - rmsd_frag2:.3f} Å")

if rmsd_complete < 2.0:
    print(f"\n🎉 Complete ligand achieves target!")