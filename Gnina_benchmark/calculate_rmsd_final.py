#!/usr/bin/env python3
# calculate_rmsd_final.py
"""
Simple and robust RMSD calculation - outputs only final result
Handles hydrogen differences automatically
"""
import sys
from rdkit import Chem
import numpy as np
import warnings

# Suppress RDKit warnings
warnings.filterwarnings('ignore')

if len(sys.argv) != 3:
    print("Uso: python calculate_rmsd_final.py <referencia.pdb> <prediccion.pdb>")
    sys.exit(1)

ref_file = sys.argv[1]
probe_file = sys.argv[2]

def read_molecule(filename):
    """Read molecule, try multiple strategies"""
    for remove_h in [False, True]:
        for sanitize in [True, False]:
            try:
                mol = Chem.MolFromPDBFile(filename, removeHs=remove_h, sanitize=sanitize)
                if mol is not None:
                    return mol
            except:
                continue
    return None

# Read molecules
ref = read_molecule(ref_file)
probe = read_molecule(probe_file)

if ref is None or probe is None:
    print("Error: Could not read structures")
    sys.exit(1)

# Always work with heavy atoms only to avoid H issues
ref_noH = Chem.RemoveHs(ref, sanitize=False)
probe_noH = Chem.RemoveHs(probe, sanitize=False)

# Check atom counts
if ref_noH.GetNumAtoms() != probe_noH.GetNumAtoms():
    print(f"Error: Atom count mismatch (Ref: {ref_noH.GetNumAtoms()}, Docked: {probe_noH.GetNumAtoms()})")
    sys.exit(1)

try:
    # Get conformers
    ref_conf = ref_noH.GetConformer()
    probe_conf = probe_noH.GetConformer()
    n = ref_noH.GetNumAtoms()
    
    # Get coordinates
    ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in range(n)])
    probe_coords = np.array([probe_conf.GetAtomPosition(i) for i in range(n)])
    
    # Center structures
    ref_center = ref_coords - ref_coords.mean(axis=0)
    probe_center = probe_coords - probe_coords.mean(axis=0)
    
    # Kabsch alignment
    H = probe_center.T @ ref_center
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Correct for reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Apply rotation
    probe_aligned = (R @ probe_center.T).T
    
    # Calculate RMSD
    diff = ref_center - probe_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    # Output in expected format
    print(f"RMSD: {rmsd:.3f} A")
    sys.exit(0)
    
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)