#!/usr/bin/env python3
"""
Test Smina with relaxed parameters - Updated
"""

import subprocess
from pathlib import Path
import numpy as np

WORK_DIR = Path(r'G:\Smina_benchmark')

def kabsch_rmsd(P, Q):
    """Calculate RMSD using Kabsch algorithm"""
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
    """Extract coordinates from PDB/PDBQT"""
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

def run_smina(seed, exhaustiveness, box_size=1.0):
    """Run Smina with specified parameters"""
    
    receptor = WORK_DIR / "receptor.pdbqt"
    ligand = WORK_DIR / "ligand.pdbqt"
    output = WORK_DIR / f"output_seed{seed}.pdbqt"
    
    # Fragment 1 coordinates
    center_x, center_y, center_z = 16.48, 15.24, 25.54
    size = 20 * box_size
    
    cmd = [
        "smina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size),
        "--size_y", str(size),
        "--size_z", str(size),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "9",
        "--energy_range", "3",
        "--seed", str(seed),
        "--out", str(output)
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
            cwd=str(WORK_DIR)
        )
        
        if output.exists() and output.stat().st_size > 0:
            return output
        else:
            return None
            
    except Exception as e:
        return None

def test_parameters():
    """Test different exhaustiveness values"""
    
    print("="*70)
    print("TESTING DIFFERENT EXHAUSTIVENESS VALUES")
    print("="*70)
    
    reference = WORK_DIR / "ligand_A1IDX_prepared.pdb"
    if not reference.exists():
        reference = WORK_DIR / "ligand.pdbqt"
    
    ref_coords = extract_coords(reference)
    if ref_coords is None:
        print("✗ Cannot load reference coordinates")
        return
    
    print(f"\nReference: {reference.name} ({len(ref_coords)} atoms)")
    
    # Test exhaustiveness values
    exh_values = [8, 16, 32, 64]
    
    for exh in exh_values:
        print(f"\n{'='*70}")
        print(f"Testing exhaustiveness = {exh}")
        print(f"{'='*70}")
        
        successes = 0
        rmsds = []
        
        # Try 5 different seeds
        for seed in range(1, 6):
            print(f"  [{seed}/5] seed={seed:2d} ", end='', flush=True)
            
            output = run_smina(seed, exh, box_size=1.0)
            
            if output and output.exists():
                docked_coords = extract_coords(output)
                
                if docked_coords is not None and len(docked_coords) >= len(ref_coords):
                    min_atoms = min(len(ref_coords), len(docked_coords))
                    rmsd = kabsch_rmsd(ref_coords[:min_atoms], docked_coords[:min_atoms])
                    
                    print(f"{rmsd:.3f}Å ✓")
                    successes += 1
                    rmsds.append(rmsd)
                else:
                    print("✗ (bad coords)")
            else:
                print("✗ (failed)")
        
        # Summary for this exhaustiveness
        print(f"\n  Results: {successes}/5 successful")
        
        if rmsds:
            print(f"  RMSD: {np.mean(rmsds):.3f} ± {np.std(rmsds):.3f} Å")
            print(f"  Best: {min(rmsds):.3f} Å")
        
        if successes >= 4:
            print(f"  ✓ Good success rate at exh={exh}")
            break
        elif successes == 0:
            print(f"  ✗ Complete failure at exh={exh}")
            print(f"  ! Check PDBQT files or box parameters")
            break

def main():
    print("="*70)
    print("Smina Parameter Testing - Updated")
    print("="*70)
    
    # Check files
    receptor = WORK_DIR / "receptor.pdbqt"
    ligand = WORK_DIR / "ligand.pdbqt"
    
    if not receptor.exists():
        print(f"✗ receptor.pdbqt not found")
        return
    
    if not ligand.exists():
        print(f"✗ ligand.pdbqt not found")
        return
    
    print(f"\n✓ receptor.pdbqt ({receptor.stat().st_size} bytes)")
    print(f"✓ ligand.pdbqt ({ligand.stat().st_size} bytes)")
    
    test_parameters()
    
    print("\n" + "="*70)
    print("RECOMMENDATION")
    print("="*70)
    print("\nUse the lowest exhaustiveness that gives good results.")
    print("Then run multiseed test with that value.")

if __name__ == "__main__":
    main()