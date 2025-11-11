#!/usr/bin/env python3
# calculate_rmsd_mcs.py
"""
RMSD calculation using MCS (Maximum Common Substructure) alignment
Solves atom ordering issues - topology-aware method
"""
import sys
from rdkit import Chem
from rdkit.Chem import rdFMCS
import numpy as np
import warnings

# Suppress RDKit warnings
warnings.filterwarnings('ignore')

if len(sys.argv) != 3:
    print("Usage: python calculate_rmsd_mcs.py <reference.pdb> <docked.pdb>")
    sys.exit(1)

ref_file = sys.argv[1]
probe_file = sys.argv[2]

def read_molecule(filename):
    """Read molecule, try multiple strategies"""
    for remove_h in [True, False]:
        for sanitize in [False, True]:
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

# Always work with heavy atoms only
ref_noH = Chem.RemoveHs(ref, sanitize=False)
probe_noH = Chem.RemoveHs(probe, sanitize=False)

try:
    # Find Maximum Common Substructure
    mcs = rdFMCS.FindMCS(
        [ref_noH, probe_noH],
        timeout=10,
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareAny
    )
    
    if mcs.numAtoms == 0:
        print("Error: No common substructure found")
        sys.exit(1)
    
    # Get atom matches based on topology
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    ref_match = ref_noH.GetSubstructMatch(mcs_mol)
    probe_match = probe_noH.GetSubstructMatch(mcs_mol)
    
    if not ref_match or not probe_match or len(ref_match) != len(probe_match):
        print("Error: Could not match substructures")
        sys.exit(1)
    
    # Get coordinates of matched atoms (topology-aware)
    ref_conf = ref_noH.GetConformer()
    probe_conf = probe_noH.GetConformer()
    
    ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in ref_match])
    probe_coords = np.array([probe_conf.GetAtomPosition(i) for i in probe_match])
    
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