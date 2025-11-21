#!/usr/bin/env python3
"""
Test Smina parameters with Docker
"""

import subprocess
from pathlib import Path
import numpy as np

WORK_DIR = Path(r'G:\Smina_benchmark')
DOCKER_IMAGE = "smina-vinardo:latest"

def kabsch_rmsd(P, Q):
    """Calculate RMSD"""
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    
    H = P_c.T @ Q_c
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    
    d = np.sign(np.linalg.det(V @ U.T))
    E = np.diag([1, 1, d])
    
    R = V @ E @ U.T
    Q_aligned = Q_c @ R.T
    
    rmsd = np.sqrt(np.mean(np.sum((P_c - Q_aligned)**2, axis=1)))
    return rmsd

def extract_coords(pdb_file):
    """Extract coordinates"""
    coords = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except:
                    continue
    return np.array(coords) if coords else None

def run_smina_docker(seed, exhaustiveness):
    """Run Smina in Docker"""
    
    output = WORK_DIR / f"output_seed{seed}.pdbqt"
    
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{WORK_DIR}:/data",
        DOCKER_IMAGE,
        "smina",
        "--receptor", "/data/receptor.pdbqt",
        "--ligand", "/data/ligand.pdbqt",
        "--center_x", "16.48",
        "--center_y", "15.24",
        "--center_z", "25.54",
        "--size_x", "20",
        "--size_y", "20",
        "--size_z", "20",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "9",
        "--energy_range", "3",
        "--seed", str(seed),
        "--out", f"/data/output_seed{seed}.pdbqt"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        return output if output.exists() and output.stat().st_size > 0 else None
    except:
        return None

def test_parameters():
    """Test different exhaustiveness"""
    
    print("="*70)
    print("TESTING WITH DOCKER")
    print("="*70)
    
    reference = WORK_DIR / "ligand_A1IDX_prepared.pdb"
    if not reference.exists():
        reference = WORK_DIR / "ligand.pdbqt"
    
    ref_coords = extract_coords(reference)
    if ref_coords is None:
        print("✗ Cannot load reference")
        return
    
    print(f"\nReference: {reference.name} ({len(ref_coords)} atoms)")
    
    exh_values = [8, 16, 32]
    
    for exh in exh_values:
        print(f"\n{'='*70}")
        print(f"Exhaustiveness = {exh}")
        print(f"{'='*70}")
        
        successes = 0
        rmsds = []
        
        for seed in range(1, 6):
            print(f"  [{seed}/5] seed={seed} ", end='', flush=True)
            
            output = run_smina_docker(seed, exh)
            
            if output:
                docked = extract_coords(output)
                
                if docked is not None and len(docked) >= len(ref_coords):
                    min_atoms = min(len(ref_coords), len(docked))
                    rmsd = kabsch_rmsd(ref_coords[:min_atoms], docked[:min_atoms])
                    
                    print(f"{rmsd:.3f}Å ✓")
                    successes += 1
                    rmsds.append(rmsd)
                else:
                    print("✗")
            else:
                print("✗")
        
        print(f"\n  Success: {successes}/5")
        
        if rmsds:
            print(f"  RMSD: {np.mean(rmsds):.3f} ± {np.std(rmsds):.3f} Å")
            print(f"  Best: {min(rmsds):.3f} Å")
        
        if successes >= 4:
            print(f"  ✓ Good at exh={exh}")
            break

def main():
    test_parameters()

if __name__ == "__main__":
    main()