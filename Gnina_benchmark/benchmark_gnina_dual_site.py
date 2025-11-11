#!/usr/bin/env python3
"""
Gnina Dual-Site Benchmark Script
Tests both binding sites (Fragment 1 and Fragment 2) separately
Optimized for multi-GPU setup (2x RTX3060) and AMD Threadripper 3690x
"""

import os
import sys
import subprocess
import numpy as np
from pathlib import Path
import time
import shutil

# ==================== CONFIGURATION ====================

WORKDIR = Path(r"G:\Gnina_benchmark")

# Docker configuration
DOCKER_IMAGE = "gnina/gnina:latest"
DOCKER_VOLUME = f"{WORKDIR}:/data"

# Two binding sites (both fragments)
SITES = {
    'Fragment_1': {
        'center_x': 16.48,
        'center_y': 15.24,
        'center_z': 25.54,
        'size_x': 20,
        'size_y': 20,
        'size_z': 20,
        'description': 'Fragment 1 binding site'
    },
    'Fragment_2': {
        'center_x': 8.44,
        'center_y': 15.48,
        'center_z': 74.27,
        'size_x': 20,
        'size_y': 20,
        'size_z': 20,
        'description': 'Fragment 2 binding site'
    }
}

# Docking parameters
DEFAULT_NUM_RUNS = 5
DEFAULT_EXHAUSTIVENESS = 32
DEFAULT_NUM_MODES = 9
DEFAULT_CPU_THREADS = 48
DEFAULT_CNN_SCORING = "rescore"

# Files
RECEPTOR = WORKDIR / "9FMM_protein.pdb"
LIGAND = WORKDIR / "ligand_A1IDX.pdb"

# ==================== FUNCTIONS ====================

def prepare_site_reference(site_name, fragment_index):
    """Prepare reference ligand for specific site"""
    from rdkit import Chem
    
    # Read original ligand with both fragments
    original = WORKDIR / "ligand_A1IDX.pdb"
    mol = Chem.MolFromPDBFile(str(original), removeHs=False, sanitize=False)
    
    if mol is None:
        print(f"   ❌ Could not read original ligand")
        return None
    
    # Get fragments
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    
    if len(frags) <= fragment_index:
        print(f"   ❌ Fragment {fragment_index} not found")
        return None
    
    # Write selected fragment as reference for this site
    ref_file = WORKDIR / f"ligand_exp_{site_name}.pdb"
    writer = Chem.PDBWriter(str(ref_file))
    writer.write(frags[fragment_index])
    writer.close()
    
    return ref_file


def run_gnina_docking(site_name, site_config, run_number):
    """Execute a single Gnina docking run for a specific site"""
    output_sdf = WORKDIR / f"output_{site_name}_{run_number}.sdf"
    log_file = WORKDIR / f"gnina_log_{site_name}_{run_number}.txt"
    
    # Build Gnina command
    cmd = [
        "docker", "run",
        "--gpus", "all",
        "-v", DOCKER_VOLUME,
        "--rm",
        DOCKER_IMAGE,
        "gnina",
        
        # Input files
        "-r", "/data/9FMM_protein.pdb",
        "-l", "/data/ligand_A1IDX.pdb",
        
        # Output
        "-o", f"/data/output_{site_name}_{run_number}.sdf",
        
        # Box definition
        "--center_x", str(site_config['center_x']),
        "--center_y", str(site_config['center_y']),
        "--center_z", str(site_config['center_z']),
        "--size_x", str(site_config['size_x']),
        "--size_y", str(site_config['size_y']),
        "--size_z", str(site_config['size_z']),
        
        # Search parameters
        "--exhaustiveness", str(DEFAULT_EXHAUSTIVENESS),
        "--num_modes", str(DEFAULT_NUM_MODES),
        "--cpu", str(DEFAULT_CPU_THREADS),
        
        # CNN scoring
        "--cnn_scoring", DEFAULT_CNN_SCORING,
        
        # Seed
        "--seed", str(run_number * 42),
    ]
    
    try:
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=3600
            )
        
        if result.returncode != 0:
            return False
        
        if not output_sdf.exists():
            return False
        
        return True
        
    except:
        return False


def convert_and_calculate_rmsd(site_name, run_number, reference_file):
    """Convert SDF to PDB and calculate RMSD"""
    from rdkit import Chem
    
    sdf_file = WORKDIR / f"output_{site_name}_{run_number}.sdf"
    pdb_file = WORKDIR / f"output_{site_name}_{run_number}.pdb"
    rmsd_file = WORKDIR / f"rmsd_{site_name}_{run_number}.txt"
    
    # Convert SDF to PDB
    try:
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
        mol = next(supplier)
        
        if mol is None:
            return None
        
        writer = Chem.PDBWriter(str(pdb_file))
        writer.write(mol)
        writer.close()
    except:
        return None
    
    # Calculate RMSD - try MCS method first (topology-aware)
    rmsd_scripts = [
        WORKDIR / "calculate_rmsd_mcs.py",     # Best: MCS topology-aware
        WORKDIR / "calculate_rmsd_final.py",   # Fallback: simple method
        WORKDIR / "calculate_rmsd.py"          # Last resort: original
    ]
    
    for rmsd_script in rmsd_scripts:
        if not rmsd_script.exists():
            continue
        
        try:
            result = subprocess.run(
                ["python", str(rmsd_script), 
                 str(reference_file), str(pdb_file)],
                capture_output=True,
                text=True,
                check=True,
                timeout=30
            )
            
            # Save and parse
            with open(rmsd_file, 'w') as f:
                f.write(result.stdout)
            
            if "RMSD:" in result.stdout:
                rmsd_value = float(result.stdout.split("RMSD:")[1].split("A")[0].strip())
                return rmsd_value
        except:
            continue  # Try next script
    
    return None


def benchmark_site(site_name, site_config, num_runs):
    """Run benchmark for a specific site"""
    print(f"\n{'='*70}")
    print(f"🎯 BENCHMARKING SITE: {site_name}")
    print(f"{'='*70}")
    print(f"Description: {site_config['description']}")
    print(f"Center: ({site_config['center_x']:.2f}, {site_config['center_y']:.2f}, {site_config['center_z']:.2f})")
    print(f"Runs: {num_runs}")
    print(f"{'='*70}\n")
    
    # Prepare reference for this site
    print(f"[PREP] Preparing reference for {site_name}...")
    fragment_idx = 0 if site_name == 'Fragment_1' else 1
    reference_file = prepare_site_reference(site_name, fragment_idx)
    
    if reference_file is None:
        print(f"   ❌ Failed to prepare reference")
        return []
    
    print(f"   ✅ Reference prepared: {reference_file.name}")
    
    # Run docking for this site
    rmsd_values = []
    site_start_time = time.time()
    
    for run in range(1, num_runs + 1):
        print(f"\n[{site_name}] Run {run}/{num_runs}")
        print("-" * 50)
        
        run_start = time.time()
        
        # Docking
        print(f"[1/3] Running Gnina docking...")
        if not run_gnina_docking(site_name, site_config, run):
            print(f"   ❌ Docking failed")
            rmsd_values.append(None)
            continue
        
        run_time = time.time() - run_start
        print(f"   ✅ Docking completed in {run_time:.1f}s")
        
        # RMSD calculation
        print(f"[2/3] Converting and calculating RMSD...")
        rmsd = convert_and_calculate_rmsd(site_name, run, reference_file)
        
        if rmsd is not None:
            print(f"   📊 RMSD: {rmsd:.3f} Å")
            rmsd_values.append(rmsd)
        else:
            print(f"   ⚠️  RMSD calculation failed")
            rmsd_values.append(None)
    
    site_time = time.time() - site_start_time
    print(f"\n⏱️  Site benchmark time: {site_time/60:.1f} minutes")
    
    return rmsd_values


def generate_site_summary(site_name, rmsd_values):
    """Generate summary for a specific site"""
    valid_rmsds = [r for r in rmsd_values if r is not None]
    
    if not valid_rmsds:
        return None
    
    summary = {
        'site': site_name,
        'successful_runs': len(valid_rmsds),
        'total_runs': len(rmsd_values),
        'average': np.mean(valid_rmsds),
        'minimum': np.min(valid_rmsds),
        'maximum': np.max(valid_rmsds),
        'std_dev': np.std(valid_rmsds),
        'success_2A': sum(1 for r in valid_rmsds if r <= 2.0),
        'success_3A': sum(1 for r in valid_rmsds if r <= 3.0),
    }
    
    return summary


def print_comparison_summary(summaries):
    """Print comparison between both sites"""
    print("\n" + "="*70)
    print("📊 DUAL-SITE BENCHMARK COMPARISON")
    print("="*70)
    
    for summary in summaries:
        if summary is None:
            continue
        
        print(f"\n🎯 {summary['site']}:")
        print(f"   Successful runs: {summary['successful_runs']}/{summary['total_runs']}")
        print(f"   Average RMSD: {summary['average']:.3f} Å")
        print(f"   Range: {summary['minimum']:.3f} - {summary['maximum']:.3f} Å")
        print(f"   Std Dev: {summary['std_dev']:.3f} Å")
        print(f"   Success ≤2Å: {summary['success_2A']}/{summary['successful_runs']} ({100*summary['success_2A']/summary['successful_runs']:.1f}%)")
        print(f"   Success ≤3Å: {summary['success_3A']}/{summary['successful_runs']} ({100*summary['success_3A']/summary['successful_runs']:.1f}%)")
    
    # Determine best site
    if all(s is not None for s in summaries) and len(summaries) == 2:
        print(f"\n🏆 COMPARISON:")
        if summaries[0]['average'] < summaries[1]['average']:
            best = summaries[0]
            diff = summaries[1]['average'] - summaries[0]['average']
        else:
            best = summaries[1]
            diff = summaries[0]['average'] - summaries[1]['average']
        
        print(f"   Best site: {best['site']}")
        print(f"   Average RMSD: {best['average']:.3f} Å")
        print(f"   Difference: {diff:.3f} Å better than other site")
        
        if best['average'] <= 2.0:
            print(f"   ✅ Excellent redocking performance!")
        elif best['average'] <= 3.0:
            print(f"   ✅ Good redocking performance")
        else:
            print(f"   ⚠️  Moderate redocking performance")
    
    print("\n" + "="*70)


def main():
    """Main dual-site benchmark execution"""
    print("\n" + "="*70)
    print("🧬 GNINA DUAL-SITE REDOCKING BENCHMARK")
    print("="*70)
    print(f"Target: 9FMM (Redocking - Both Sites)")
    print(f"Hardware: 2x RTX3060 GPUs + {DEFAULT_CPU_THREADS} CPU threads")
    
    # Get number of runs
    if len(sys.argv) > 1:
        num_runs = int(sys.argv[1])
    else:
        num_runs = DEFAULT_NUM_RUNS
    
    print(f"Runs per site: {num_runs}")
    print(f"Exhaustiveness: {DEFAULT_EXHAUSTIVENESS}")
    print(f"Total docking runs: {num_runs * len(SITES)}")
    print("="*70)
    
    # Check prerequisites
    if not RECEPTOR.exists() or not LIGAND.exists():
        print("❌ Error: Input files not found")
        sys.exit(1)
    
    # Benchmark each site
    all_summaries = []
    total_start = time.time()
    
    for site_name, site_config in SITES.items():
        rmsd_values = benchmark_site(site_name, site_config, num_runs)
        summary = generate_site_summary(site_name, rmsd_values)
        all_summaries.append(summary)
    
    # Total time
    total_time = time.time() - total_start
    print(f"\n⏱️  Total benchmark time: {total_time/60:.1f} minutes")
    
    # Print comparison
    print_comparison_summary(all_summaries)
    
    print(f"\n📁 Results saved in: {WORKDIR}")
    print(f"   Files: output_Fragment_*_*.sdf, rmsd_Fragment_*_*.txt")
    print("\n✅ Dual-site benchmark completed successfully!")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n⚠️  Benchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)