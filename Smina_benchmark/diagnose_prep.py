#!/usr/bin/env python3
"""
Diagnose Smina preparation and create missing receptor.pdbqt
"""

from pathlib import Path
import subprocess

WORK_DIR = Path(r'G:\Smina_benchmark')

def find_files():
    """Find ligand and receptor files"""
    
    print("="*70)
    print("FINDING FILES")
    print("="*70)
    
    # Ligand files (prioritize prepared versions)
    ligand_candidates = [
        "ligand_A1IDX_prepared.pdb",
        "ligand_A1IDX_fragment1_prepared.pdb",
        "ligand_A1IDX.pdb"
    ]
    
    ligand_file = None
    for candidate in ligand_candidates:
        path = WORK_DIR / candidate
        if path.exists():
            ligand_file = path
            print(f"✓ Ligand: {candidate}")
            break
    
    if not ligand_file:
        print("✗ No ligand file found")
        return None, None
    
    # Receptor file
    receptor_candidates = [
        "9FMM_protein.pdb",
        "receptor.pdb",
        "protein.pdb"
    ]
    
    receptor_file = None
    for candidate in receptor_candidates:
        # Check parent directories too
        for check_dir in [WORK_DIR, WORK_DIR.parent]:
            path = check_dir / candidate
            if path.exists():
                receptor_file = path
                print(f"✓ Receptor: {path}")
                break
        if receptor_file:
            break
    
    if not receptor_file:
        print("✗ No receptor file found")
    
    return ligand_file, receptor_file

def diagnose_pdb(pdb_file, is_ligand=True):
    """Diagnose PDB file"""
    
    file_type = "LIGAND" if is_ligand else "RECEPTOR"
    
    print(f"\n{'='*70}")
    print(f"DIAGNOSING {file_type}: {pdb_file.name}")
    print("="*70)
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # Count atoms
    atom_count = 0
    hetatm_count = 0
    hydrogen_count = 0
    elements = set()
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_count += 1
            atom_name = line[12:16].strip()
            if atom_name.startswith('H'):
                hydrogen_count += 1
        elif line.startswith('HETATM'):
            hetatm_count += 1
        
        if line.startswith('ATOM') or line.startswith('HETATM'):
            # Extract element (columns 77-78)
            if len(line) > 77:
                element = line[76:78].strip()
                if element:
                    elements.add(element)
    
    print(f"\nAtom statistics:")
    print(f"  ATOM records: {atom_count}")
    print(f"  HETATM records: {hetatm_count}")
    print(f"  Hydrogens: {hydrogen_count}")
    print(f"  Elements: {sorted(elements)}")
    
    # Check for issues
    issues = []
    
    if is_ligand:
        if hetatm_count == 0 and atom_count == 0:
            issues.append("No atoms found")
        if hydrogen_count == 0:
            issues.append("No hydrogens (may need addition)")
    else:
        if atom_count == 0:
            issues.append("No ATOM records")
        if hydrogen_count == 0:
            issues.append("No hydrogens (recommend adding polar H)")
    
    if issues:
        print(f"\n! Issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print(f"\n✓ Structure looks OK")
    
    return len(issues) == 0

def convert_to_pdbqt(pdb_file, output_file, is_ligand=True):
    """Convert PDB to PDBQT using OpenBabel"""
    
    print(f"\n{'='*70}")
    print(f"CONVERTING: {pdb_file.name} → {output_file.name}")
    print("="*70)
    
    obabel_path = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
    
    if not Path(obabel_path).exists():
        print(f"✗ OpenBabel not found at: {obabel_path}")
        return False
    
    if is_ligand:
        # Ligand: add hydrogens, compute charges, flexible
        cmd = [
            obabel_path,
            str(pdb_file),
            "-O", str(output_file),
            "-h",  # Add hydrogens
            "--partialcharge", "gasteiger",  # Gasteiger charges
            "-p", "7.4"  # pH
        ]
    else:
        # Receptor: rigid, polar hydrogens only
        cmd = [
            obabel_path,
            str(pdb_file),
            "-O", str(output_file),
            "-xr",  # Rigid
            "-p", "7.4"  # Add polar H at pH 7.4
        ]
    
    print(f"\nCommand: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if output_file.exists():
            size = output_file.stat().st_size
            print(f"✓ Created: {output_file.name} ({size} bytes)")
            
            # Quick validation
            with open(output_file, 'r') as f:
                content = f.read()
            
            atom_count = content.count('ATOM')
            hetatm_count = content.count('HETATM')
            
            print(f"  Atoms: {atom_count + hetatm_count}")
            
            if is_ligand:
                has_root = 'ROOT' in content
                has_torsdof = 'TORSDOF' in content
                
                print(f"  ROOT: {'✓' if has_root else '✗ MISSING'}")
                print(f"  TORSDOF: {'✓' if has_torsdof else '✗ MISSING'}")
                
                if not has_root or not has_torsdof:
                    print(f"\n! WARNING: PDBQT incomplete for docking")
                    print(f"! Try using MGLTools prepare_ligand4.py instead")
                    return False
            
            return True
        else:
            print(f"✗ Output file not created")
            return False
            
    except subprocess.CalledProcessError as e:
        print(f"✗ Conversion failed:")
        print(f"  {e.stderr}")
        return False

def create_test_command():
    """Create test Smina command with relaxed parameters"""
    
    print(f"\n{'='*70}")
    print("TEST SMINA COMMAND")
    print("="*70)
    
    cmd = """
# Test with RELAXED parameters first:

smina \\
  --receptor receptor.pdbqt \\
  --ligand ligand.pdbqt \\
  --center_x 16.48 --center_y 15.24 --center_z 25.54 \\
  --size_x 20 --size_y 20 --size_z 20 \\
  --exhaustiveness 8 \\
  --num_modes 9 \\
  --energy_range 3 \\
  --seed 1 \\
  --out result.pdbqt \\
  --log result.log

# If that works, gradually increase exhaustiveness:
# 8 → 16 → 32 → 64 → 128
"""
    
    print(cmd)

def main():
    # Find files
    ligand_pdb, receptor_pdb = find_files()
    
    if not ligand_pdb:
        print("\n✗ Cannot proceed without ligand file")
        return
    
    # Diagnose
    ligand_ok = diagnose_pdb(ligand_pdb, is_ligand=True)
    
    if receptor_pdb:
        receptor_ok = diagnose_pdb(receptor_pdb, is_ligand=False)
    else:
        print("\n✗ Receptor PDB not found")
        receptor_ok = False
    
    # Convert to PDBQT
    print(f"\n{'='*70}")
    print("CONVERSION TO PDBQT")
    print("="*70)
    
    ligand_pdbqt = WORK_DIR / "ligand.pdbqt"
    receptor_pdbqt = WORK_DIR / "receptor.pdbqt"
    
    if ligand_pdb:
        print(f"\nConverting ligand...")
        convert_to_pdbqt(ligand_pdb, ligand_pdbqt, is_ligand=True)
    
    if receptor_pdb:
        print(f"\nConverting receptor...")
        convert_to_pdbqt(receptor_pdb, receptor_pdbqt, is_ligand=False)
    
    # Create test command
    create_test_command()
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY & NEXT STEPS")
    print("="*70)
    
    print(f"\nFiles created:")
    if ligand_pdbqt.exists():
        print(f"  ✓ {ligand_pdbqt}")
    else:
        print(f"  ✗ ligand.pdbqt")
    
    if receptor_pdbqt.exists():
        print(f"  ✓ {receptor_pdbqt}")
    else:
        print(f"  ✗ receptor.pdbqt")
    
    print(f"\nTest with:")
    print(f"  1. Run the smina command above")
    print(f"  2. If it works, increase exhaustiveness")
    print(f"  3. Previous failure at exh=256 was TOO HIGH")
    print(f"  4. Start with exh=8, then 16, 32, 64")

if __name__ == "__main__":
    main()