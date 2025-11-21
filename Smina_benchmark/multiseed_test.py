#!/usr/bin/env python3
"""
Multi-seed runner: 10 runs con mismo config, diferentes seeds
"""
import subprocess
import shutil
from pathlib import Path
import numpy as np
import time

class MultiSeedRunner:
    def __init__(self):
        self.base = Path(r"G:\Smina_benchmark")
        self.results_dir = self.base / "multiseed"
        self.results_dir.mkdir(exist_ok=True)
        
        self.protein = self.base / "9FMM_protein.pdb"
        self.ligand = self.base / "ligand_A1IDX_fragment2.pdb"
        self.ref_coords = self._extract(self.ligand)
    
    def _extract(self, pdb):
        coords = []
        with open(pdb, 'r') as f:
            for line in f:
                if line.startswith('ENDMDL'):
                    break
                if line.startswith(('ATOM', 'HETATM')) and not line[12:16].strip().startswith('H'):
                    try:
                        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                    except:
                        pass
        return np.array(coords)
    
    def _rmsd(self, P, Q):
        P_c = P - P.mean(axis=0)
        Q_c = Q - Q.mean(axis=0)
        C = np.dot(Q_c.T, P_c)
        V, S, W = np.linalg.svd(C)
        if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
            V[:, -1] = -V[:, -1]
        Q_rot = np.dot(Q_c - Q_c.mean(axis=0), np.dot(V, W))
        return np.sqrt(np.mean(np.sum((P_c - Q_rot)**2, axis=1)))
    
    def run_seed(self, box, exh, modes, seed, run_id):
        run_dir = self.results_dir / f"run_{run_id:03d}"
        run_dir.mkdir(exist_ok=True)
        
        shutil.copy(self.protein, run_dir / "protein.pdb")
        shutil.copy(self.ligand, run_dir / "ligand.pdb")
        
        print(f"[{run_id:03d}] seed={seed:<10} ", end="", flush=True)
        
        cmd = [
            "docker", "run", "--rm",
            "-v", f"{run_dir}:/workspace", "-w", "/workspace",
            "smina-vinardo:latest", "smina",
            "-r", "protein.pdb", "-l", "ligand.pdb",
            "--autobox_ligand", "ligand.pdb",
            "--autobox_add", str(box),
            "--exhaustiveness", str(exh),
            "--num_modes", str(modes),
            "--seed", str(seed),
            "--scoring", "vinardo",
            "-o", "docked.pdb"
        ]
        
        try:
            subprocess.run(cmd, capture_output=True, timeout=600)
            docked = self._extract(run_dir / "docked.pdb")
            rmsd = self._rmsd(self.ref_coords, docked)
            print(f"{rmsd:.3f}Å")
            return rmsd
        except:
            print("✗")
            return 999.0
    
    def run(self):
        # Mejor config conocida
        box, exh, modes = 1.0, 256, 50
        
        print(f"Running 20 times with different seeds")
        print(f"Config: box={box}, exh={exh}, modes={modes}\n")
        
        seeds = list(range(1, 21))  # Seeds 1-20
        rmsds = []
        
        for i, seed in enumerate(seeds):
            rmsd = self.run_seed(box, exh, modes, seed, i)
            if rmsd < 999:
                rmsds.append(rmsd)
        
        if rmsds:
            print(f"\n{'='*60}")
            print(f"STATISTICS (n={len(rmsds)}):")
            print(f"{'='*60}")
            print(f"Best:   {min(rmsds):.3f} Å")
            print(f"Worst:  {max(rmsds):.3f} Å")
            print(f"Mean:   {np.mean(rmsds):.3f} Å")
            print(f"Median: {np.median(rmsds):.3f} Å")
            print(f"Std:    {np.std(rmsds):.3f} Å")
            
            if min(rmsds) < 2.0:
                print(f"\n🎉 Found seed with RMSD < 2.0 Å!")

if __name__ == "__main__":
    runner = MultiSeedRunner()
    runner.run()