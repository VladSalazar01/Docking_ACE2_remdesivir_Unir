#!/usr/bin/env python3
"""
Fix ligand.pdbqt - regenerate with proper format for Smina
"""

from pathlib import Path
import subprocess

WORK_DIR = Path(r'G:\Smina_benchmark')

def inspect_ligand():
    """Show the problematic lines"""
    
    ligand = WORK_DIR / "ligand.pdbqt"
    
    print("="*70)
    print("INSPECTING LIGAND.PDBQT")
    print("="*70)
    
    with open(ligand, 'r') as f:
        lines = f.readlines()
    
    print(f"\nTotal lines: {len(lines)}")
    print(f"\nLines 105-120 (around error line 112):")
    
    for i in range(104, min(120, len(lines))):
        marker = " <-- LINE 112" if i == 111 else ""
        print(f"  {i+1:3d}: {lines[i].rstrip()}{marker}")
    
    print(f"\nAll lines:")
    for i, line in enumerate(lines):
        print(f"  {i+1:3d}: {line.rstrip()}")

def regenerate_ligand_simple():
    """Regenerate ligand with minimal OpenBabel flags"""
    
    print(f"\n{'='*70}")
    print("REGENERATING LIGAND.PDBQT")
    print("="*70)
    
    # Use prepared PDB as input
    input_pdb = WORK_DIR / "ligand_A1IDX_prepared.pdb"
    if not input_pdb.exists():
        input_pdb = WORK_DIR / "ligand_A1IDX.pdb"
    
    if not input_pdb.exists():
        print(f"✗ No input PDB found")
        return False
    
    output_pdbqt = WORK_DIR / "ligand_fixed.pdbqt"
    
    obabel = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"
    
    # Try minimal conversion - just format change
    print(f"\n1. Trying minimal conversion (no hydrogen addition)...")
    
    cmd = [
        obabel,
        str(input_pdb),
        "-O", str(output_pdbqt),
        "-xn"  # No bond typing
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if output_pdbqt.exists():
            print(f"✓ Created {output_pdbqt.name}")
            
            # Check structure
            with open(output_pdbqt, 'r') as f:
                content = f.read()
            
            has_root = 'ROOT' in content
            has_torsdof = 'TORSDOF' in content
            
            print(f"  ROOT: {has_root}")
            print(f"  TORSDOF: {has_torsdof}")
            
            if not has_root or not has_torsdof:
                print(f"\n  ! Missing ROOT/TORSDOF")
                print(f"  ! OpenBabel doesn't add these")
                print(f"  ! Need MGLTools or manual addition")
                return False
            
            return True
        
    except Exception as e:
        print(f"✗ Failed: {e}")
        return False

def add_pdbqt_header():
    """Manually add ROOT/TORSDOF to make valid PDBQT"""
    
    print(f"\n{'='*70}")
    print("ADDING PDBQT HEADER MANUALLY")
    print("="*70)
    
    # Read current ligand
    ligand = WORK_DIR / "ligand.pdbqt"
    
    with open(ligand, 'r') as f:
        lines = f.readlines()
    
    # Remove any problematic lines
    clean_lines = []
    for i, line in enumerate(lines):
        # Skip lines with unknown tags
        if any(tag in line for tag in ['MODEL', 'ENDMDL', 'CONECT', 'MASTER', 'END']):
            print(f"  Removing line {i+1}: {line.rstrip()}")
            continue
        clean_lines.append(line)
    
    # Build proper PDBQT
    output_lines = []
    
    # Add REMARK
    output_lines.append("REMARK  Name = ligand\n")
    output_lines.append("REMARK  0 active torsions:\n")
    output_lines.append("REMARK  status: ('A' for Active; 'I' for Inactive)\n")
    
    # Add ROOT
    output_lines.append("ROOT\n")
    
    # Add ATOM lines
    for line in clean_lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            output_lines.append(line)
    
    # Add ENDROOT
    output_lines.append("ENDROOT\n")
    
    # Add TORSDOF (0 for rigid ligand to test)
    output_lines.append("TORSDOF 0\n")
    
    # Write fixed version
    output = WORK_DIR / "ligand_manual.pdbqt"
    
    with open(output, 'w') as f:
        f.writelines(output_lines)
    
    print(f"\n✓ Created: {output.name}")
    print(f"  Lines: {len(output_lines)}")
    
    print(f"\nFirst 20 lines:")
    for i, line in enumerate(output_lines[:20]):
        print(f"  {i+1:2d}: {line.rstrip()}")
    
    return output

def test_fixed_ligand(ligand_file):
    """Test Smina with fixed ligand"""
    
    print(f"\n{'='*70}")
    print("TESTING FIXED LIGAND")
    print("="*70)
    
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{WORK_DIR}:/data",
        "smina-vinardo:latest",
        "smina",
        "--receptor", "/data/receptor.pdbqt",
        "--ligand", f"/data/{ligand_file.name}",
        "--center_x", "16.48",
        "--center_y", "15.24",
        "--center_z", "25.54",
        "--size_x", "20",
        "--size_y", "20",
        "--size_z", "20",
        "--exhaustiveness", "8",
        "--num_modes", "9",
        "--seed", "1",
        "--out", "/data/test_fixed.pdbqt"
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        print(f"\nSTDOUT:")
        print(result.stdout)
        
        if result.stderr:
            print(f"\nSTDERR:")
            print(result.stderr)
        
        output = WORK_DIR / "test_fixed.pdbqt"
        if output.exists() and output.stat().st_size > 0:
            print(f"\n✓ SUCCESS! Docking worked with fixed ligand")
            return True
        else:
            print(f"\n✗ Still failed")
            return False
            
    except Exception as e:
        print(f"\n✗ Error: {e}")
        return False

def main():
    inspect_ligand()
    
    # Try manual fix first
    fixed_ligand = add_pdbqt_header()
    
    if fixed_ligand:
        success = test_fixed_ligand(fixed_ligand)
        
        if success:
            print(f"\n{'='*70}")
            print("SOLUTION")
            print("="*70)
            print(f"\nUse: {fixed_ligand.name}")
            print(f"Update scripts to use this file instead of ligand.pdbqt")
        else:
            print(f"\n{'='*70}")
            print("ALTERNATIVE SOLUTION")
            print("="*70)
            print(f"\nTry using Vina instead of Smina")
            print(f"Or use Gnina (already working)")

if __name__ == "__main__":
    main()