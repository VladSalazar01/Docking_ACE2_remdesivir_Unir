#!/usr/bin/env python3
"""
DIAGNÓSTICO: ¿POR QUÉ LIGANDO COMPLETO FALLA?
=============================================
Investiga por qué ligand_A1IDX.pdb (56 átomos) da peor RMSD
que ligand_A1IDX_fragment2.pdb (28 átomos)
"""

import subprocess
from pathlib import Path
from typing import Dict, List, Tuple
import re

class LigandDiagnostics:
    """Diagnóstico de ligandos"""
    
    def __init__(self, base_dir: str = r"G:\Smina_benchmark"):
        self.base_dir = Path(base_dir)
        self.fragment2 = self.base_dir / "ligand_A1IDX_fragment2.pdb"
        self.complete = self.base_dir / "ligand_A1IDX.pdb"
        self.protein = self.base_dir / "9FMM_protein.pdb"
        
    def analyze_pdb_structure(self, pdb_file: Path) -> Dict:
        """Analiza estructura de un PDB"""
        print(f"\n{'=' * 70}")
        print(f"ANALYZING: {pdb_file.name}")
        print(f"{'=' * 70}")
        
        if not pdb_file.exists():
            print(f"✗ File not found: {pdb_file}")
            return {}
        
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
        
        # Parse atoms
        atoms = []
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_data = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'resName': line[17:20].strip(),
                    'chainID': line[21:22].strip(),
                    'resSeq': int(line[22:26].strip()),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'element': line[76:78].strip()
                }
                atoms.append(atom_data)
        
        # Statistics
        total_atoms = len(atoms)
        elements = [a['element'] for a in atoms if a['element']]
        
        from collections import Counter
        element_counts = Counter(elements)
        
        # Heavy atoms (not H)
        heavy_atoms = [a for a in atoms if a['element'] != 'H']
        hydrogens = [a for a in atoms if a['element'] == 'H']
        
        # Center of mass
        com_x = sum(a['x'] for a in heavy_atoms) / len(heavy_atoms)
        com_y = sum(a['y'] for a in heavy_atoms) / len(heavy_atoms)
        com_z = sum(a['z'] for a in heavy_atoms) / len(heavy_atoms)
        
        # Bounding box
        min_x = min(a['x'] for a in atoms)
        max_x = max(a['x'] for a in atoms)
        min_y = min(a['y'] for a in atoms)
        max_y = max(a['y'] for a in atoms)
        min_z = min(a['z'] for a in atoms)
        max_z = max(a['z'] for a in atoms)
        
        size_x = max_x - min_x
        size_y = max_y - min_y
        size_z = max_z - min_z
        
        info = {
            'total_atoms': total_atoms,
            'heavy_atoms': len(heavy_atoms),
            'hydrogens': len(hydrogens),
            'elements': dict(element_counts),
            'center_of_mass': (com_x, com_y, com_z),
            'bounding_box': ((min_x, max_x), (min_y, max_y), (min_z, max_z)),
            'size': (size_x, size_y, size_z)
        }
        
        # Print results
        print(f"\n📊 STRUCTURE INFO:")
        print(f"  Total atoms: {total_atoms}")
        print(f"  Heavy atoms: {len(heavy_atoms)}")
        print(f"  Hydrogens: {len(hydrogens)}")
        print(f"\n🧪 COMPOSITION:")
        for elem, count in element_counts.most_common():
            print(f"  {elem}: {count}")
        print(f"\n📍 SPATIAL INFO:")
        print(f"  Center of mass: ({com_x:.2f}, {com_y:.2f}, {com_z:.2f})")
        print(f"  Size (x, y, z): ({size_x:.2f}, {size_y:.2f}, {size_z:.2f}) Å")
        print(f"  Volume (approx): {size_x * size_y * size_z:.2f} Å³")
        
        return info
    
    def compare_ligands(self):
        """Compara fragment2 vs complete"""
        print(f"\n{'#' * 70}")
        print("LIGAND COMPARISON: Fragment2 vs Complete")
        print(f"{'#' * 70}")
        
        frag2_info = self.analyze_pdb_structure(self.fragment2)
        complete_info = self.analyze_pdb_structure(self.complete)
        
        if not frag2_info or not complete_info:
            return
        
        print(f"\n{'=' * 70}")
        print("COMPARISON SUMMARY")
        print(f"{'=' * 70}")
        
        # Difference in atoms
        diff_total = complete_info['total_atoms'] - frag2_info['total_atoms']
        diff_heavy = complete_info['heavy_atoms'] - frag2_info['heavy_atoms']
        
        print(f"\n📈 SIZE DIFFERENCE:")
        print(f"  Total atoms: {complete_info['total_atoms']} - {frag2_info['total_atoms']} = +{diff_total}")
        print(f"  Heavy atoms: {complete_info['heavy_atoms']} - {frag2_info['heavy_atoms']} = +{diff_heavy}")
        
        # Size difference
        frag2_vol = frag2_info['size'][0] * frag2_info['size'][1] * frag2_info['size'][2]
        complete_vol = complete_info['size'][0] * complete_info['size'][1] * complete_info['size'][2]
        vol_ratio = complete_vol / frag2_vol if frag2_vol > 0 else 0
        
        print(f"\n📦 VOLUME DIFFERENCE:")
        print(f"  Fragment2: {frag2_vol:.2f} Å³")
        print(f"  Complete: {complete_vol:.2f} Å³")
        print(f"  Ratio: {vol_ratio:.2f}x larger")
        
        # Distance between centers
        com_frag2 = frag2_info['center_of_mass']
        com_complete = complete_info['center_of_mass']
        
        dist = ((com_complete[0] - com_frag2[0])**2 + 
                (com_complete[1] - com_frag2[1])**2 + 
                (com_complete[2] - com_frag2[2])**2)**0.5
        
        print(f"\n📍 CENTER OF MASS SHIFT:")
        print(f"  Distance: {dist:.2f} Å")
        
        if dist > 2.0:
            print(f"  ⚠️  WARNING: Large shift! Complete ligand center is very different")
        
        print(f"\n{'=' * 70}")
        print("🔍 INTERPRETATION")
        print(f"{'=' * 70}")
        
        # Interpretation
        if diff_heavy > 15:
            print(f"\n1. ⚠️  LARGE SIZE DIFFERENCE:")
            print(f"   Complete has {diff_heavy} more heavy atoms than Fragment2")
            print(f"   This could cause:")
            print(f"   - Steric clashes with protein")
            print(f"   - Difficulty fitting in binding pocket")
            print(f"   - Need for induced fit (flexible receptor)")
        
        if vol_ratio > 2.0:
            print(f"\n2. ⚠️  VOLUME MISMATCH:")
            print(f"   Complete is {vol_ratio:.1f}x larger by volume")
            print(f"   With box_add=1 Å, Complete might not fit!")
        
        if dist > 2.0:
            print(f"\n3. ⚠️  CENTER MISMATCH:")
            print(f"   Centers differ by {dist:.2f} Å")
            print(f"   AutoBox might be centering incorrectly")
        
        print(f"\n{'=' * 70}")
        print("💡 RECOMMENDATIONS")
        print(f"{'=' * 70}")
        
        print(f"\n✓ Fragment2 is likely the CORRECT ligand because:")
        print(f"  1. It gives better RMSD (3.2 Å vs 5.8 Å)")
        print(f"  2. Smaller molecules are easier to dock accurately")
        print(f"  3. It probably matches the crystallized ligand")
        
        print(f"\n✓ To improve Complete ligand docking:")
        print(f"  1. Use larger box (autobox_add=2-3)")
        print(f"  2. Try flexible receptor")
        print(f"  3. Verify it's the correct chemical structure")
        print(f"  4. Check if Fragment2 IS the complete ligand")
        
        print(f"\n✓ Verify against PDB:")
        print(f"  1. Download 9FMM.cif from RCSB PDB")
        print(f"  2. Extract the crystallized ligand")
        print(f"  3. Compare atom count with Fragment2")
    
    def test_box_sizes_complete(self):
        """Prueba diferentes box sizes para ligando completo"""
        print(f"\n{'=' * 70}")
        print("TESTING BOX SIZES FOR COMPLETE LIGAND")
        print(f"{'=' * 70}")
        
        if not self.complete.exists():
            print(f"✗ Complete ligand not found: {self.complete}")
            return
        
        box_sizes = [1.5, 2.0, 2.5, 3.0]
        results = []
        
        for box_add in box_sizes:
            print(f"\n▶ Testing box_add = {box_add} Å...")
            
            output_file = self.base_dir / f"test_complete_box{box_add}.pdb"
            log_file = self.base_dir / f"test_complete_box{box_add}.log"
            
            cmd = [
                "smina",
                "-r", str(self.protein),
                "-l", str(self.complete),
                "--autobox_ligand", str(self.complete),
                "--autobox_add", str(box_add),
                "--exhaustiveness", "64",
                "--num_modes", "20",
                "--scoring", "vinardo",
                "-o", str(output_file),
                "--log", str(log_file)
            ]
            
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=600
                )
                
                if result.returncode == 0:
                    # Parse RMSD
                    with open(log_file, 'r') as f:
                        for line in f:
                            if line.strip().startswith('1 '):
                                parts = line.split()
                                if len(parts) >= 3:
                                    rmsd = float(parts[2])
                                    print(f"  RMSD: {rmsd:.3f} Å")
                                    results.append((box_add, rmsd))
                                    break
                else:
                    print(f"  ✗ Failed: {result.stderr[:100]}")
                    
            except subprocess.TimeoutExpired:
                print(f"  ✗ Timeout")
            except Exception as e:
                print(f"  ✗ Error: {e}")
        
        if results:
            print(f"\n{'=' * 70}")
            print("RESULTS SUMMARY")
            print(f"{'=' * 70}")
            
            best_box, best_rmsd = min(results, key=lambda x: x[1])
            
            for box_add, rmsd in sorted(results):
                marker = " ← BEST" if box_add == best_box else ""
                print(f"  Box {box_add:.1f} Å: RMSD {rmsd:.3f} Å{marker}")
            
            print(f"\n💡 CONCLUSION:")
            if best_rmsd < frag2_info.get('best_rmsd', 3.187):
                print(f"  ✓ Complete CAN match Fragment2 with box={best_box}")
            else:
                print(f"  ✗ Complete still worse than Fragment2 (3.187 Å)")
                print(f"  → Fragment2 is probably the correct ligand")


def main():
    """Main execution"""
    print("""
╔═══════════════════════════════════════════════════════════════════════════════╗
║                    LIGAND DIAGNOSTICS: Fragment2 vs Complete                  ║
║                                                                               ║
║  Goal: Understand why Complete (56 atoms) gives worse RMSD                  ║
║        than Fragment2 (28 atoms)                                             ║
╚═══════════════════════════════════════════════════════════════════════════════╝
    """)
    
    diag = LigandDiagnostics()
    
    # Comparison
    diag.compare_ligands()
    
    # Test different box sizes for complete
    print(f"\n\n{'#' * 70}")
    print("OPTIONAL: Test Complete with larger boxes?")
    print(f"{'#' * 70}")
    response = input("Run box size test? (y/n): ").lower()
    
    if response == 'y':
        diag.test_box_sizes_complete()
    else:
        print("\nSkipped. You can run this test later by:")
        print("  python diagnose_complete_ligand.py")
    
    print("\n✅ Diagnostics complete!")


if __name__ == "__main__":
    main()
