#!/usr/bin/env python3
"""
Gnina Benchmark Script for Redocking - Enhanced Version with Config Support
Optimized for multi-GPU setup (2x RTX3060) and AMD Threadripper 3690x
Author: Automated Pipeline
Date: 2025-11-06
Version: 2.0 (with config file support)
"""

import os
import sys
import subprocess
import glob
import numpy as np
from pathlib import Path
import time
import configparser

# ==================== DEFAULT CONFIGURATION ====================

# Directories and files
DEFAULT_WORKDIR = Path(r"G:\Gnina_benchmark")
DEFAULT_RECEPTOR = "9FMM_protein.pdb"
DEFAULT_LIGAND = "ligand_A1IDX.pdb"
DEFAULT_LIGAND_EXP = "ligand_exp.pdb"

# Docker configuration
DEFAULT_DOCKER_IMAGE = "gnina/gnina:latest"

# Docking parameters - centered on Fragment 1 (active site)
DEFAULT_CENTER_X = 16.48
DEFAULT_CENTER_Y = 15.24
DEFAULT_CENTER_Z = 25.54
DEFAULT_SIZE_X = 20
DEFAULT_SIZE_Y = 20
DEFAULT_SIZE_Z = 20

# Performance parameters
DEFAULT_NUM_RUNS = 5
DEFAULT_EXHAUSTIVENESS = 32
DEFAULT_NUM_MODES = 9
DEFAULT_CPU_THREADS = 48
DEFAULT_GPUS = "0,1"
DEFAULT_CNN_SCORING = "rescore"
# Note: Gnina uses its built-in CNN model automatically
DEFAULT_TIMEOUT = 3600
DEFAULT_SEED_BASE = 42

# ==================== CONFIGURATION LOADER ====================

def load_config():
    """Load configuration from gnina_config.ini if it exists"""
    config_file = DEFAULT_WORKDIR / "gnina_config.ini"
    
    # Default configuration
    config = {
        'workdir': DEFAULT_WORKDIR,
        'receptor': DEFAULT_RECEPTOR,
        'ligand': DEFAULT_LIGAND,
        'ligand_exp': DEFAULT_LIGAND_EXP,
        'docker_image': DEFAULT_DOCKER_IMAGE,
        'center_x': DEFAULT_CENTER_X,
        'center_y': DEFAULT_CENTER_Y,
        'center_z': DEFAULT_CENTER_Z,
        'size_x': DEFAULT_SIZE_X,
        'size_y': DEFAULT_SIZE_Y,
        'size_z': DEFAULT_SIZE_Z,
        'num_runs': DEFAULT_NUM_RUNS,
        'exhaustiveness': DEFAULT_EXHAUSTIVENESS,
        'num_modes': DEFAULT_NUM_MODES,
        'cpu_threads': DEFAULT_CPU_THREADS,
        'gpus': DEFAULT_GPUS,
        'cnn_scoring': DEFAULT_CNN_SCORING,
        'timeout': DEFAULT_TIMEOUT,
        'seed_base': DEFAULT_SEED_BASE,
    }
    
    # Load from config file if it exists
    if config_file.exists():
        print(f"📄 Loading configuration from: {config_file}")
        parser = configparser.ConfigParser()
        parser.read(config_file)
        
        # Override defaults with config file values
        if 'DOCKING' in parser:
            config['num_runs'] = parser.getint('DOCKING', 'num_runs', fallback=DEFAULT_NUM_RUNS)
            config['exhaustiveness'] = parser.getint('DOCKING', 'exhaustiveness', fallback=DEFAULT_EXHAUSTIVENESS)
            config['num_modes'] = parser.getint('DOCKING', 'num_modes', fallback=DEFAULT_NUM_MODES)
            config['cpu_threads'] = parser.getint('DOCKING', 'cpu_threads', fallback=DEFAULT_CPU_THREADS)
            config['gpus'] = parser.get('DOCKING', 'gpus', fallback=DEFAULT_GPUS)
        
        if 'CNN_SCORING' in parser:
            config['cnn_scoring'] = parser.get('CNN_SCORING', 'cnn_scoring', fallback=DEFAULT_CNN_SCORING)
        
        if 'BOX' in parser:
            config['center_x'] = parser.getfloat('BOX', 'center_x', fallback=DEFAULT_CENTER_X)
            config['center_y'] = parser.getfloat('BOX', 'center_y', fallback=DEFAULT_CENTER_Y)
            config['center_z'] = parser.getfloat('BOX', 'center_z', fallback=DEFAULT_CENTER_Z)
            config['size_x'] = parser.getfloat('BOX', 'size_x', fallback=DEFAULT_SIZE_X)
            config['size_y'] = parser.getfloat('BOX', 'size_y', fallback=DEFAULT_SIZE_Y)
            config['size_z'] = parser.getfloat('BOX', 'size_z', fallback=DEFAULT_SIZE_Z)
        
        if 'PATHS' in parser:
            workdir_str = parser.get('PATHS', 'workdir', fallback=str(DEFAULT_WORKDIR))
            config['workdir'] = Path(workdir_str)
            config['receptor'] = parser.get('PATHS', 'receptor', fallback=DEFAULT_RECEPTOR)
            config['ligand'] = parser.get('PATHS', 'ligand', fallback=DEFAULT_LIGAND)
            config['ligand_exp'] = parser.get('PATHS', 'ligand_exp', fallback=DEFAULT_LIGAND_EXP)
        
        if 'DOCKER' in parser:
            config['docker_image'] = parser.get('DOCKER', 'image', fallback=DEFAULT_DOCKER_IMAGE)
        
        if 'ADVANCED' in parser:
            config['timeout'] = parser.getint('ADVANCED', 'timeout', fallback=DEFAULT_TIMEOUT)
            config['seed_base'] = parser.getint('ADVANCED', 'seed_base', fallback=DEFAULT_SEED_BASE)
        
        print("   ✅ Configuration loaded from file")
    else:
        print("ℹ️  Using default configuration (no gnina_config.ini found)")
    
    return config

# ==================== MAIN SCRIPT ====================

def check_prerequisites(config):
    """Verify all required files and tools exist"""
    print("\n🔍 Checking prerequisites...")
    
    WORKDIR = config['workdir']
    RECEPTOR = WORKDIR / config['receptor']
    LIGAND = WORKDIR / config['ligand']
    
    # Check working directory
    if not WORKDIR.exists():
        print(f"❌ Error: Working directory not found: {WORKDIR}")
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
            ["docker", "images", "-q", config['docker_image']],
            capture_output=True, text=True, check=True
        )
        if not result.stdout.strip():
            print(f"   📥 Pulling Docker image: {config['docker_image']}")
            subprocess.run(["docker", "pull", config['docker_image']], check=True)
        print(f"   ✅ Docker image ready: {config['docker_image']}")
    except subprocess.CalledProcessError:
        print(f"❌ Error: Could not verify/pull Docker image")
        sys.exit(1)
    
    # Prepare experimental ligand pose if it doesn't exist
    LIGAND_EXP = WORKDIR / config['ligand_exp']
    if not LIGAND_EXP.exists():
        print(f"   📋 Creating experimental pose: {LIGAND_EXP}")
        import shutil
        shutil.copy(LIGAND, LIGAND_EXP)
    
    print("✅ All prerequisites met!\n")


def run_gnina_docking(run_number, config):
    """Execute a single Gnina docking run"""
    WORKDIR = config['workdir']
    output_sdf = WORKDIR / f"output_{run_number}.sdf"
    log_file = WORKDIR / f"gnina_log_{run_number}.txt"
    
    # Build Gnina command
    cmd = [
        "docker", "run",
        "--gpus", "all",
        "-v", f"{WORKDIR}:/data",
        "--rm",
        config['docker_image'],
        "gnina",  # Explicit gnina executable
        
        # Input files
        "-r", f"/data/{config['receptor']}",
        "-l", f"/data/{config['ligand']}",
        
        # Output
        "-o", f"/data/output_{run_number}.sdf",
        
        # Box definition
        "--center_x", str(config['center_x']),
        "--center_y", str(config['center_y']),
        "--center_z", str(config['center_z']),
        "--size_x", str(config['size_x']),
        "--size_y", str(config['size_y']),
        "--size_z", str(config['size_z']),
        
        # Search parameters
        "--exhaustiveness", str(config['exhaustiveness']),
        "--num_modes", str(config['num_modes']),
        "--cpu", str(config['cpu_threads']),
        
        # CNN scoring
        "--cnn_scoring", config['cnn_scoring'],
        
        # Additional options
        "--seed", str(run_number * config['seed_base']),
    ]
    
    print(f"   🚀 Running docking (GPU: {config['gpus']}, CPU: {config['cpu_threads']} threads)...")
    
    try:
        with open(log_file, 'w') as log:
            result = subprocess.run(
                cmd,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=config['timeout']
            )
        
        if result.returncode != 0:
            print(f"   ⚠️  Warning: Docking may have issues (exit code: {result.returncode})")
            return False
        
        if not output_sdf.exists():
            print(f"   ❌ Error: Output file not created")
            return False
        
        return True
        
    except subprocess.TimeoutExpired:
        print(f"   ❌ Error: Docking timeout (>{config['timeout']}s)")
        return False
    except Exception as e:
        print(f"   ❌ Error: {e}")
        return False


def convert_sdf_to_pdb(run_number, config):
    """Convert SDF output to PDB for RMSD calculation"""
    WORKDIR = config['workdir']
    sdf_file = WORKDIR / f"output_{run_number}.sdf"
    pdb_file = WORKDIR / f"output_{run_number}.pdb"
    
    try:
        from rdkit import Chem
        
        supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False, sanitize=False)
        mol = next(supplier)
        
        if mol is None:
            print(f"   ⚠️  Warning: Could not read SDF file")
            return False
        
        # Sanitize if possible
        try:
            Chem.SanitizeMol(mol)
        except:
            pass  # Continue even if sanitization fails
        
        writer = Chem.PDBWriter(str(pdb_file))
        writer.write(mol)
        writer.close()
        
        return True
        
    except Exception as e:
        print(f"   ⚠️  Warning: SDF to PDB conversion failed: {e}")
        
        try:
            subprocess.run(
                ["obabel", "-isdf", str(sdf_file), "-opdb", 
                 "-O", str(pdb_file), "-d"],
                capture_output=True,
                check=True
            )
            return True
        except:
            return False


def calculate_rmsd(run_number, config):
    """Calculate RMSD using the provided calculate_rmsd.py script"""
    WORKDIR = config['workdir']
    pdb_file = WORKDIR / f"output_{run_number}.pdb"
    rmsd_file = WORKDIR / f"rmsd_{run_number}.txt"
    ligand_exp = WORKDIR / config['ligand_exp']
    
    if not pdb_file.exists():
        print(f"   ⚠️  Warning: PDB file not found, skipping RMSD")
        return None
    
    # Try RMSD scripts in order of preference
    rmsd_scripts = [
        WORKDIR / "calculate_rmsd_final.py",     # Simplest, most robust
        WORKDIR / "calculate_rmsd_smart.py",     # Smart with MCS
        WORKDIR / "calculate_rmsd_improved.py",  # Improved version
        WORKDIR / "calculate_rmsd.py"            # Original
    ]
    
    for script in rmsd_scripts:
        if not script.exists():
            continue
            
        try:
            result = subprocess.run(
                ["python", str(script), 
                 str(ligand_exp), str(pdb_file)],
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
                return None
                
        except subprocess.CalledProcessError as e:
            # Try next script
            if script == rmsd_scripts[-1]:  # Last script
                print(f"   ⚠️  Warning: RMSD calculation failed")
                if e.stderr:
                    print(f"      Error: {e.stderr}")
                return None
            continue
        except Exception as e:
            print(f"   ⚠️  Warning: RMSD calculation failed: {e}")
            return None
    
    print(f"   ⚠️  Warning: No RMSD calculation script found")
    return None


def generate_summary(rmsd_values, config):
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
    
    success_2A = sum(1 for r in valid_rmsds if r <= 2.0)
    success_3A = sum(1 for r in valid_rmsds if r <= 3.0)
    
    print(f"\n🎯 Success Rate:")
    print(f"   • RMSD ≤ 2.0 Å: {success_2A}/{len(valid_rmsds)} ({100*success_2A/len(valid_rmsds):.1f}%)")
    print(f"   • RMSD ≤ 3.0 Å: {success_3A}/{len(valid_rmsds)} ({100*success_3A/len(valid_rmsds):.1f}%)")
    
    print("\n" + "="*60)
    print(f"📁 Results saved in: {config['workdir']}")
    print("="*60 + "\n")


def main():
    """Main benchmark execution"""
    # Load configuration
    config = load_config()
    
    # Allow command-line override of num_runs
    if len(sys.argv) > 1:
        config['num_runs'] = int(sys.argv[1])
        print(f"📝 Command line override: {config['num_runs']} runs")
    
    print("\n" + "="*60)
    print("🧬 GNINA REDOCKING BENCHMARK")
    print("="*60)
    print(f"Target: 9FMM (Redocking)")
    print(f"Hardware: 2x RTX3060 GPUs + {config['cpu_threads']} CPU threads")
    print(f"Runs: {config['num_runs']}")
    print(f"Exhaustiveness: {config['exhaustiveness']}")
    print(f"CNN Scoring: {config['cnn_scoring']}")
    print("="*60 + "\n")
    
    # Check prerequisites
    check_prerequisites(config)
    
    # Store RMSD values
    rmsd_values = []
    
    # Benchmark loop
    start_time = time.time()
    
    for run in range(1, config['num_runs'] + 1):
        print(f"\n{'='*60}")
        print(f"🔄 Run {run}/{config['num_runs']}")
        print(f"{'='*60}")
        
        run_start = time.time()
        
        print(f"[1/3] Running Gnina docking...")
        if not run_gnina_docking(run, config):
            print(f"   ❌ Docking failed for run {run}")
            rmsd_values.append(None)
            continue
        
        run_time = time.time() - run_start
        print(f"   ✅ Docking completed in {run_time:.1f}s")
        
        print(f"[2/3] Converting SDF to PDB...")
        if not convert_sdf_to_pdb(run, config):
            print(f"   ❌ Conversion failed for run {run}")
            rmsd_values.append(None)
            continue
        print(f"   ✅ Conversion successful")
        
        print(f"[3/3] Calculating RMSD...")
        rmsd = calculate_rmsd(run, config)
        rmsd_values.append(rmsd)
    
    # Generate summary
    total_time = time.time() - start_time
    print(f"\n⏱️  Total benchmark time: {total_time/60:.1f} minutes")
    
    generate_summary(rmsd_values, config)
    
    print("✅ Benchmark completed successfully!")


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