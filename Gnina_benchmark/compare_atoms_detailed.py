#!/usr/bin/env python3
"""
Detailed atom-by-atom comparison to find the 3 missing/extra atoms
"""

from pathlib import Path
from rdkit import Chem
from collections import Counter

WORKDIR = Path(r"G:\Gnina_benchmark")

print("="*70)
print("🔬 DETAILED ATOM COMPARISON")
print("="*70)

# Read reference
ref_file = WORKDIR / "ligand_exp.pdb"
print(f"\n[1/3] Reference: {ref_file.name}")
ref_mol = Chem.MolFromPDBFile(str(ref_file), removeHs=False, sanitize=False)

if ref_mol:
    ref_atoms_all = [atom.GetSymbol() for atom in ref_mol.GetAtoms()]
    ref_composition_all = Counter(ref_atoms_all)
    
    print(f"   Total atoms: {len(ref_atoms_all)}")
    print(f"   Composition: {dict(ref_composition_all)}")
    
    # Heavy atoms only
    ref_mol_heavy = Chem.RemoveHs(ref_mol, sanitize=False)
    ref_atoms_heavy = [atom.GetSymbol() for atom in ref_mol_heavy.GetAtoms()]
    ref_composition_heavy = Counter(ref_atoms_heavy)
    
    print(f"   Heavy atoms: {len(ref_atoms_heavy)}")
    print(f"   Composition: {dict(ref_composition_heavy)}")
else:
    print("   ❌ Could not read file")
    exit(1)

# Read docked
docked_file = WORKDIR / "output_1.pdb"
print(f"\n[2/3] Docked: {docked_file.name}")
docked_mol = Chem.MolFromPDBFile(str(docked_file), removeHs=False, sanitize=False)

if docked_mol:
    docked_atoms_all = [atom.GetSymbol() for atom in docked_mol.GetAtoms()]
    docked_composition_all = Counter(docked_atoms_all)
    
    print(f"   Total atoms: {len(docked_atoms_all)}")
    print(f"   Composition: {dict(docked_composition_all)}")
    
    # Heavy atoms only
    docked_mol_heavy = Chem.RemoveHs(docked_mol, sanitize=False)
    docked_atoms_heavy = [atom.GetSymbol() for atom in docked_mol_heavy.GetAtoms()]
    docked_composition_heavy = Counter(docked_atoms_heavy)
    
    print(f"   Heavy atoms: {len(docked_atoms_heavy)}")
    print(f"   Composition: {dict(docked_composition_heavy)}")
else:
    print("   ❌ Could not read file")
    exit(1)

# Compare
print(f"\n[3/3] Comparison:")
print(f"\n   ALL ATOMS:")
print(f"   Reference: {ref_composition_all}")
print(f"   Docked:    {docked_composition_all}")

# Find differences
all_elements = set(ref_composition_all.keys()) | set(docked_composition_all.keys())
print(f"\n   Differences:")
has_diff = False
for element in sorted(all_elements):
    ref_count = ref_composition_all.get(element, 0)
    dock_count = docked_composition_all.get(element, 0)
    if ref_count != dock_count:
        diff = dock_count - ref_count
        symbol = "+" if diff > 0 else ""
        print(f"      {element}: Ref={ref_count}, Docked={dock_count} ({symbol}{diff})")
        has_diff = True

if not has_diff:
    print("      No differences in composition!")

# Heavy atoms comparison
print(f"\n   HEAVY ATOMS ONLY:")
print(f"   Reference: {ref_composition_heavy}")
print(f"   Docked:    {docked_composition_heavy}")

heavy_elements = set(ref_composition_heavy.keys()) | set(docked_composition_heavy.keys())
print(f"\n   Differences:")
has_heavy_diff = False
for element in sorted(heavy_elements):
    ref_count = ref_composition_heavy.get(element, 0)
    dock_count = docked_composition_heavy.get(element, 0)
    if ref_count != dock_count:
        diff = dock_count - ref_count
        symbol = "+" if diff > 0 else ""
        print(f"      {element}: Ref={ref_count}, Docked={dock_count} ({symbol}{diff})")
        has_heavy_diff = True

if not has_heavy_diff:
    print("      ✅ Heavy atoms match perfectly!")

# Check SDF directly
print(f"\n[4/4] Checking original SDF file...")
sdf_file = WORKDIR / "output_1.sdf"
if sdf_file.exists():
    supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
    sdf_mol = next(supplier, None)
    
    if sdf_mol:
        sdf_atoms_all = [atom.GetSymbol() for atom in sdf_mol.GetAtoms()]
        sdf_composition_all = Counter(sdf_atoms_all)
        
        print(f"   SDF total atoms: {len(sdf_atoms_all)}")
        print(f"   SDF composition: {dict(sdf_composition_all)}")
        
        # Heavy atoms
        sdf_mol_heavy = Chem.RemoveHs(sdf_mol, sanitize=False)
        sdf_atoms_heavy = [atom.GetSymbol() for atom in sdf_mol_heavy.GetAtoms()]
        sdf_composition_heavy = Counter(sdf_atoms_heavy)
        
        print(f"   SDF heavy atoms: {len(sdf_atoms_heavy)}")
        print(f"   SDF heavy composition: {dict(sdf_composition_heavy)}")
        
        # Compare SDF heavy with reference heavy
        print(f"\n   SDF vs Reference (heavy atoms):")
        if sdf_composition_heavy == ref_composition_heavy:
            print(f"      ✅ Compositions MATCH!")
        else:
            print(f"      ❌ Compositions DIFFER:")
            all_heavy = set(sdf_composition_heavy.keys()) | set(ref_composition_heavy.keys())
            for element in sorted(all_heavy):
                sdf_count = sdf_composition_heavy.get(element, 0)
                ref_count = ref_composition_heavy.get(element, 0)
                if sdf_count != ref_count:
                    diff = sdf_count - ref_count
                    symbol = "+" if diff > 0 else ""
                    print(f"         {element}: SDF={sdf_count}, Ref={ref_count} ({symbol}{diff})")

print("\n" + "="*70)
print("📊 DIAGNOSIS")
print("="*70)

if has_heavy_diff:
    print("\n❌ ISSUE: Heavy atom composition differs between reference and docked")
    print("   This means Gnina changed the molecular structure")
    print("   Possible causes:")
    print("   - Gnina added/removed atoms during docking")
    print("   - Input ligand has issues")
    print("   - Protonation state changes")
else:
    print("\n✅ Heavy atoms match! The difference is only in hydrogens")
    print("   Solution: Use removeHs=True in RMSD calculation")

print("="*70)