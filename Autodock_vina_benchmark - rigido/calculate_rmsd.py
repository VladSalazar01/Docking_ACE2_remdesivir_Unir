# calculate_rmsd.py
import sys
from rdkit import Chem
import numpy as np

if len(sys.argv) != 3:
    print("Uso: python calculate_rmsd.py <referencia.pdb> <prediccion.pdb>")
    sys.exit(1)

ref_file = sys.argv[1]
probe_file = sys.argv[2]

try:
    ref = Chem.MolFromPDBFile(ref_file, removeHs=False)
    probe = Chem.MolFromPDBFile(probe_file, removeHs=False)
    
    if ref is None:
        print(f"Error: No se pudo leer {ref_file}")
        sys.exit(1)
    
    if probe is None:
        print(f"Error: No se pudo leer {probe_file}")
        sys.exit(1)
    
    if ref.GetNumAtoms() != probe.GetNumAtoms():
        print(f"Error: Ref={ref.GetNumAtoms()} vs Probe={probe.GetNumAtoms()}")
        sys.exit(1)
    
    # Obtener coordenadas
    ref_conf = ref.GetConformer()
    probe_conf = probe.GetConformer()
    n = ref.GetNumAtoms()
    
    ref_coords = np.array([ref_conf.GetAtomPosition(i) for i in range(n)])
    probe_coords = np.array([probe_conf.GetAtomPosition(i) for i in range(n)])
    
    # Centrar estructuras
    ref_center = ref_coords - ref_coords.mean(axis=0)
    probe_center = probe_coords - probe_coords.mean(axis=0)
    
    # Alineamiento Kabsch
    H = probe_center.T @ ref_center
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    
    # Corregir reflexión
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    
    # Aplicar rotación
    probe_aligned = (R @ probe_center.T).T
    
    # Calcular RMSD
    diff = ref_center - probe_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    print(f"RMSD: {rmsd:.3f} A")
    
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)