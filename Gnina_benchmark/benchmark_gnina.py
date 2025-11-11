#!/usr/bin/env python3
"""
Gnina Benchmark Script for Redocking
Optimized for multi-GPU setup (2x RTX3060) and AMD Threadripper 3690x
Author: Automated Pipeline
Date: 2025-11-06
"""

import os
import sys
import subprocess
import glob
import numpy as np
from pathlib import Path
import time

# ==================== CONFIGURATION ====================

# Directories and files
WORKDIR = Path(r"G:\Gnina_benchmark")
RECEPTOR = WORKDIR / "9FMM_protein.pdb"
LIGAND = WORKDIR / "ligand_A1IDX.pdb"
LIGAND_EXP = WORKDIR / "ligand_exp.pdb"
CALCULATE_RMSD = WORKDIR / "calculate_rmsd.py"

# Docker configuration
DOCKER_IMAGE = "gnina/gnina:latest"
DOCKER_VOLUME = f"{WORKDIR}:/data"

# Docking box parameters - centered on Fragment 1 (active site)
CENTER_X = 16.48
CENTER_Y = 15.24
CENTER_Z = 25.54
SIZE_X = 20
SIZE_Y = 20
SIZE_Z = 20

# Gnina parameters optimized for your hardware
NUM_RUNS = 5  # Default number of benchmark runs (easily adjustable)
EXHAUSTIVENESS = 32  # Search effort (increase for better results, slower)
NUM_MODES = 9  # Number of binding modes to generate
CPU_THREADS = 48  # Your CPU thread count
GPUS = "0,1"  # Both RTX3060 GPUs

# CNN scoring (Gnina's key feature for better accuracy)
CNN_SCORING = "rescore"  # Options: "rescore", "refinement", "all", "none"
# Note: Gnina uses its built-in CNN model automatically, no need to specify --cnn

# ==================== FUNCTIONS ====================

def check_prerequisites():
    """Verify all required files and tools exist"""
    print("🔍 Checking prerequisites...")
    
    # Check working directory
    if not WORKDIR.exists():
        print(f"❌ Error: Working directory not found: {WORKDIR}")
        print(f"   Please create the directory and place input files there.")
        sys.exit(1)
    
    # Check input files
    required_files = [RECEPTOR, LIGAND]
    for file in required_files:
        if not file.exists():
            print(f"❌ Error: Required file not found: {file}")
            sys.exit(1)
    
    # Check Docker
    try:
        subprocess.run(["docker", "--version"], 
                      capture_output=True, check=True)
        print("   ✅ Docker found")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("❌ Error: Docker not found or not running")
        sys.exit(1)
    
    # Check if Docker image exists or pull it
    try:
        result = subprocess.run(
            ["docker", "images", "-q", DOCKER_IMAGE],
            capture_output=True, text=True, check=True
        )
        if not result.stdout.strip():
            print(f"   📥 Pulling Docker image: {DOCKER_IMAGE}")
            subprocess.run(["docker", "pull", DOCKER_IMAGE], check=True)
        print(f"   ✅ Docker image ready: {DOCKER_IMAGE}")
    except subprocess.CalledProcessError:
        print(f"❌ Error: Could not verify/pull Docker image")
        sys.exit(1)
    
    # Prepare experimental ligand pose if it doesn't exist
    if not LIGAND_EXP.exists():
        print(f"   📋 Creating experimental pose: {LIGAND_EXP}")
        import shutil
        shutil.copy(LIGAND, LIGAND_EXP)
    
    print("✅ All prerequisites met!\n")


def run_gnina_docking(run_number):
    """Execute a single Gnina docking run"""
    output_sdf = WORKDIR / f"output_{run_number}.sdf"
    log_file = WORKDIR / f"gnina_log_{run_number}.txt"
    
    # Build Gnina command
    cmd = [
        "docker", "run",
        "--gpus", "all",  # Enable GPU access
        "-v", DOCKER_VOLUME,
        "--rm",
        DOCKER_IMAGE,
        "gnina",  # Explicit gnina executable
        
        # Input files
        "-r", "/data/9FMM_protein.pdb",
        "-l", "/data/ligand_A1IDX.pdb",
        
        # Output
        "-o", f"/data/output_{run_number}.sdf",
        
        # Box definition
        "--center_x", str(CENTER_X),
        "--center_y", str(CENTER_Y),
        "--center_z", str(CENTER_Z),
        "--size_x", str(SIZE_X),
        "--size_y", str(SIZE_Y),
        "--size_z", str(SIZE_Z),
        
        # Search parameters
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--cpu", str(CPU_THREADS),
        
        # CNN scoring (Gnina's advantage over Vina)
        "--cnn_scoring", CNN_SCORING,
        
        # Additional options
        "--seed", str(run_number * 42),  # Different seed per run
    ]
    
    print(f"   🚀 Running docking (GPU: {GPUS}, CPU: {CPU_THREADS} threads)...")
    
    try:
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=3600  # 1 hour timeout
            )
        
        if result.returncode != 0:
            print(f"   ⚠️  Warning: Docking may have issues (exit code: {result.returncode})")
            return False
        
        if not output_sdf.exists():
            print(f"   ❌ Error: Output file not created")
            return False
        
        return True
        
    except subprocess.TimeoutExpired:
        print(f"   ❌ Error: Docking timeout (>1 hour)")
        return False
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return False


def convert_sdf_to_pdb(run_number):
    """Convert SDF output to PDB for RMSD calculation"""
    sdf_file = WORKDIR / f"output_{run_number}.sdf"
    pdb_file = WORKDIR / f"output_{run_number}.pdb"
    
    try:
        # Use RDKit to convert
        from rdkit import Chem
        
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
        mol = next(supplier)  # Get first pose (best scored)
        
        if mol is None:
            print(f"   ⚠️  Warning: Could not read SDF file")
            return False
        
        # Sanitize if possible
        try:
            Chem.SanitizeMol(mol)
        except:
            pass  # Continue even if sanitization fails
        
        # Write PDB
        writer = Chem.PDBWriter(str(pdb_file))
        writer.write(mol)
        writer.close()
        
        return True
        
    except Exception as e:
        print(f"   ⚠️  Warning: SDF to PDB conversion failed: {e}")
        
        # Fallback to obabel if RDKit fails
        try:
            subprocess.run(
                ["obabel", "-isdf", str(sdf_file), "-opdb", 
                 "-O", str(pdb_file), "-d"],  # -d removes hydrogens
                capture_output=True,
                check=True
            )
            return True
        except:
            return False


def calculate_rmsd(run_number):
    """Calculate RMSD using the provided calculate_rmsd.py script"""
    pdb_file = WORKDIR / f"output_{run_number}.pdb"
    rmsd_file = WORKDIR / f"rmsd_{run_number}.txt"
    
    if not pdb_file.exists():
        print(f"   ⚠️  Warning: PDB file not found, skipping RMSD")
        return None
    
    # Try RMSD scripts in order of preference
    rmsd_scripts = [
        WORKDIR / "calculate_rmsd_mcs.py",       # Best: Topology-aware MCS
        WORKDIR / "calculate_rmsd_final.py",     # Good: Simple heavy atom
        WORKDIR / "calculate_rmsd_smart.py",     # Backup: Smart with MCS
        WORKDIR / "calculate_rmsd_improved.py",  # Backup: Improved version
        WORKDIR / "calculate_rmsd.py"            # Last: Original
    ]
    
    for script in rmsd_scripts:
        if not script.exists():
            continue
            
        try:
            result = subprocess.run(
                ["python", str(script), 
                 str(LIGAND_EXP), str(pdb_file)],
                capture_output=True,
                text=True,
                check=True
            )
            
            # Save RMSD output
            with open(rmsd_file, 'w') as f:
                f.write(result.stdout)
            
            # Parse RMSD value
            if "RMSD:" in result.stdout:
                rmsd_value = float(result.stdout.split("RMSD:")[1].split("A")[0].strip())
                print(f"   📊 {result.stdout.strip()}")
                return rmsd_value
            else:
                print(f"   ⚠️  Warning: Could not parse RMSD from output")
                print(f"      Output: {result.stdout[:200]}")
                return None
                
        except subprocess.CalledProcessError as e:
            # Show detailed error info
            if script == rmsd_scripts[-1]:  # Last script
                print(f"   ⚠️  Warning: All RMSD scripts failed")
                print(f"      Last error from {script.name}:")
                if e.stdout:
                    print(f"      STDOUT: {e.stdout[:200]}")
                if e.stderr:
                    print(f"      STDERR: {e.stderr[:200]}")
                return None
            continue
        except Exception as e:
            print(f"   ⚠️  Warning: RMSD calculation failed: {e}")
            return None
    
    print(f"   ⚠️  Warning: No RMSD calculation script found")
    return None


def generate_summary(rmsd_values):
    """Generate statistical summary of benchmark results"""
    print("\n" + "="*60)
    print("📈 BENCHMARK SUMMARY")
    print("="*60)
    
    valid_rmsds = [r for r in rmsd_values if r is not None]
    
    if not valid_rmsds:
        print("❌ No valid RMSD values calculated")
        return
    
    print(f"\n✅ Successful runs: {len(valid_rmsds)}/{len(rmsd_values)}")
    print(f"\n📊 RMSD Statistics:")
    print(f"   • Average: {np.mean(valid_rmsds):.3f} Å")
    print(f"   • Minimum: {np.min(valid_rmsds):.3f} Å")
    print(f"   • Maximum: {np.max(valid_rmsds):.3f} Å")
    print(f"   • Std Dev: {np.std(valid_rmsds):.3f} Å")
    
    # Success rate based on RMSD threshold
    success_2A = sum(1 for r in valid_rmsds if r <= 2.0)
    success_3A = sum(1 for r in valid_rmsds if r <= 3.0)
    
    print(f"\n🎯 Success Rate:")
    print(f"   • RMSD ≤ 2.0 Å: {success_2A}/{len(valid_rmsds)} ({100*success_2A/len(valid_rmsds):.1f}%)")
    print(f"   • RMSD ≤ 3.0 Å: {success_3A}/{len(valid_rmsds)} ({100*success_3A/len(valid_rmsds):.1f}%)")
    
    print("\n" + "="*60)
    print(f"📁 Results saved in: {WORKDIR}")
    print("="*60 + "\n")


def main():
    """Main benchmark execution"""
    print("\n" + "="*60)
    print("🧬 GNINA REDOCKING BENCHMARK")
    print("="*60)
    print(f"Target: 9FMM (Redocking)")
    print(f"Hardware: 2x RTX3060 GPUs + {CPU_THREADS} CPU threads")
    print(f"Runs: {NUM_RUNS}")
    print(f"Exhaustiveness: {EXHAUSTIVENESS}")
    print("="*60 + "\n")
    
    # Check prerequisites
    check_prerequisites()
    
    # Store RMSD values for summary
    rmsd_values = []
    
    # Benchmark loop
    start_time = time.time()
    
    for run in range(1, NUM_RUNS + 1):
        print(f"\n{'='*60}")
        print(f"🔄 Run {run}/{NUM_RUNS}")
        print(f"{'='*60}")
        
        run_start = time.time()
        
        # Step 1: Docking
        print(f"[1/3] Running Gnina docking...")
        if not run_gnina_docking(run):
            print(f"   ❌ Docking failed for run {run}")
            rmsd_values.append(None)
            continue
        
        run_time = time.time() - run_start
        print(f"   ✅ Docking completed in {run_time:.1f}s")
        
        # Step 2: Convert to PDB
        print(f"[2/3] Converting SDF to PDB...")
        if not convert_sdf_to_pdb(run):
            print(f"   ❌ Conversion failed for run {run}")
            rmsd_values.append(None)
            continue
        print(f"   ✅ Conversion successful")
        
        # Step 3: Calculate RMSD
        print(f"[3/3] Calculating RMSD...")
        rmsd = calculate_rmsd(run)
        rmsd_values.append(rmsd)
    
    # Generate summary
    total_time = time.time() - start_time
    print(f"\n⏱️  Total benchmark time: {total_time/60:.1f} minutes")
    
    generate_summary(rmsd_values)
    
    print("✅ Benchmark completed successfully!")


if __name__ == "__main__":
    try:
        # Allow custom number of runs via command line
        if len(sys.argv) > 1:
            NUM_RUNS = int(sys.argv[1])
            print(f"📝 Custom run count: {NUM_RUNS}")
        
        main()
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Benchmark interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ Fatal error: {e}")
        sys.exit(1)