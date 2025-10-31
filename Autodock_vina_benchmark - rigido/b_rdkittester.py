from rdkit import Chem
from rdkit.Chem import AllChem

# Verificar archivos
ref = Chem.MolFromPDBFile('ligand_exp.pdb', removeHs=True)
print(f"Referencia cargada: {ref is not None}")
if ref:
    print(f"  Átomos: {ref.GetNumAtoms()}")

probe = Chem.MolFromPDBFile('output_1.pdb', removeHs=True)
print(f"Pose predicha cargada: {probe is not None}")
if probe:
    print(f"  Átomos: {probe.GetNumAtoms()}")

# Si ambos se cargan, calcular RMSD
if ref and probe:
    ref = Chem.RemoveHs(ref)
    probe = Chem.RemoveHs(probe)
    rmsd = AllChem.GetBestRMS(ref, probe)
    print(f"RMSD: {rmsd:.3f} Å")