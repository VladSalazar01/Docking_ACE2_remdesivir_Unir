from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

ref = Chem.MolFromPDBFile('ligand_exp_clean.pdb', removeHs=False)
probe = Chem.MolFromPDBFile('output_1_test.pdb', removeHs=False)

print(f"Referencia: {ref.GetNumAtoms()} átomos")
print(f"Predicción: {probe.GetNumAtoms()} átomos")

if ref.GetNumAtoms() == probe.GetNumAtoms():
    # Intentar alineamiento tradicional primero
    try:
        rmsd = AllChem.GetBestRMS(ref, probe)
        print(f"\n✓ RMSD (con alineamiento): {rmsd:.3f} Å")
    except RuntimeError as e:
        print(f"\n⚠ GetBestRMS falló: {e}")
        print("Usando RMSD sin alineamiento (coordenadas directas)...")
        
        # Obtener coordenadas
        ref_conf = ref.GetConformer()
        probe_conf = probe.GetConformer()
        
        ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in range(ref.GetNumAtoms())])
        probe_coords = np.array([probe_conf.GetAtomPosition(i) for i in range(probe.GetNumAtoms())])
        
        # RMSD simple (sin alineamiento)
        rmsd_simple = np.sqrt(np.mean(np.sum((ref_coords - probe_coords)**2, axis=1)))
        print(f"✓ RMSD (sin alineamiento): {rmsd_simple:.3f} Å")
        
        # Intentar alineamiento manual con Kabsch
        # Centrar ambas estructuras
        ref_centered = ref_coords - ref_coords.mean(axis=0)
        probe_centered = probe_coords - probe_coords.mean(axis=0)
        
        # Matriz de covarianza
        H = probe_centered.T @ ref_centered
        
        # SVD
        U, S, Vt = np.linalg.svd(H)
        
        # Matriz de rotación
        R = Vt.T @ U.T
        
        # Asegurar rotación correcta (det = 1)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Aplicar rotación
        probe_aligned = (R @ probe_centered.T).T
        
        # RMSD después de alineamiento óptimo
        rmsd_aligned = np.sqrt(np.mean(np.sum((ref_centered - probe_aligned)**2, axis=1)))
        print(f"✓ RMSD (con alineamiento Kabsch): {rmsd_aligned:.3f} Å")
else:
    print(f"\n⚠ Diferentes números de átomos: {ref.GetNumAtoms()} vs {probe.GetNumAtoms()}")