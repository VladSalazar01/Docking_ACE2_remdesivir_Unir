#!/usr/bin/env python3
"""
Reproducir EXACTAMENTE la configuración que dio 3.187 Å
"""
import subprocess
import shutil
from pathlib import Path
import numpy as np

class RMSDCalculator:
    @staticmethod
    def extract_coords(pdb_file):
        coords = []
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ENDMDL'):
                    break
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    if line[12:16].strip().startswith('H'):
                        continue
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
                    except:
                        continue
        return np.array(coords)
    
    @staticmethod
    def kabsch_rmsd(P, Q):
        if len(P) != len(Q):
            return 999.0
        P_c = P - P.mean(axis=0)
        Q_c = Q - Q.mean(axis=0)
        C = np.dot(Q_c.T, P_c)
        V, S, W = np.linalg.svd(C)
        if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
            V[:, -1] = -V[:, -1]
        U = np.dot(V, W)
        Q_rot = np.dot(Q_c, U)
        return np.sqrt(np.mean(np.sum((P_c - Q_rot)**2, axis=1)))

def test_exact_config():
    """
    Según tu descripción:
    - ligand_A1IDX_fragment2.pdb
    - AUTOBOX_ADD=1
    - EXHAUSTIVENESS=128
    - SCORING=vinardo
    """
    
    base = Path(r"G:\Smina_benchmark")
    test_dir = base / "test_exact_318"
    test_dir.mkdir(exist_ok=True)
    
    protein = base / "9FMM_protein.pdb"
    ligand = base / "ligand_A1IDX_fragment2.pdb"
    
    shutil.copy(protein, test_dir / "protein.pdb")
    shutil.copy(ligand, test_dir / "ligand.pdb")
    
    print("Testing EXACT configuration from your 3.187 Å result:")
    print("  - Fragment2")
    print("  - Autobox add: 1")
    print("  - Exhaustiveness: 128")
    print("  - Scoring: vinardo")
    print()
    
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{test_dir}:/workspace",
        "-w", "/workspace",
        "smina-vinardo:latest", "smina",
        "-r", "protein.pdb",
        "-l", "ligand.pdb",
        "--autobox_ligand", "ligand.pdb",
        "--autobox_add", "1",
        "--exhaustiveness", "128",
        "--num_modes", "40",
        "--scoring", "vinardo",
        "-o", "docked.pdb",
        "--log", "smina.log"
    ]
    
    print("Running...")
    result = subprocess.run(docker_cmd, capture_output=True, text=True, timeout=300)
    
    if result.returncode != 0:
        print(f"✗ Failed: {result.stderr}")
        return
    
    print("✓ Completed\n")
    
    # Calcular RMSD con Kabsch
    ref_coords = RMSDCalculator.extract_coords(ligand)
    docked_coords = RMSDCalculator.extract_coords(test_dir / "docked.pdb")
    
    rmsd = RMSDCalculator.kabsch_rmsd(ref_coords, docked_coords)
    
    print(f"{'='*60}")
    print(f"RESULT:")
    print(f"{'='*60}")
    print(f"RMSD: {rmsd:.3f} Å")
    print(f"Expected: 3.187 Å")
    print(f"Difference: {abs(rmsd - 3.187):.3f} Å")
    
    if abs(rmsd - 3.187) < 0.5:
        print("\n✓ REPRODUCED! (within 0.5 Å)")
    else:
        print(f"\n⚠ Different result")
        print(f"\nPossible reasons:")
        print(f"  1. Reference ligand is different")
        print(f"  2. Protein structure changed")
        print(f"  3. Smina version difference")
        print(f"  4. Random seed variation")
    
    # Mostrar log
    print(f"\n{'='*60}")
    print("SMINA LOG (first mode):")
    print(f"{'='*60}")
    with open(test_dir / "smina.log") as f:
        in_results = False
        for line in f:
            if 'mode |' in line:
                in_results = True
            if in_results:
                print(line.rstrip())
                if line.startswith('3 '):
                    break

if __name__ == "__main__":
    test_exact_config()