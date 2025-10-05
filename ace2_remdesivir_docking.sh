#!/bin/bash
# Perform molecular docking of remdesivir to ACE2 (PDB 6M0J) using GNINA
# Requirements: Miniconda with 'docking' environment (Python 3.9, RDKit, PyMOL-open-source, Open Babel), Docker with gnina/gnina:latest, NVIDIA drivers/CUDA
# Outputs: Receptor (6M0J_clean.pdbqt), ligand (remdesivir_aligned.pdbqt), docking results (test_gnina_specific.pdbqt)

# Define variables
set PDB_ID="6M0J"
set CHAIN="A"  # ACE2 chain in 6M0J
set SMILES="CC[C@@H]1C(=O)N(C(=O)N(C2=C1C(=N)NC(N)=C2C#N)CC(O)CO)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O"
set ACTIVE_SITE_RESIDUES="30+31+34+83+353"

# Download PDB file from RCSB
wget https://files.rcsb.org/download/${PDB_ID}.pdb

# Verify chain and residues in original PDB
grep "^ATOM" ${PDB_ID}.pdb | grep " ${CHAIN} " | head -n 20 > chain_${CHAIN}.txt
grep " ASP ${CHAIN}  30" ${PDB_ID}.pdb >> residues.txt || echo "Asp30 not found" >> residues.txt
grep " LYS ${CHAIN}  31" ${PDB_ID}.pdb >> residues.txt || echo "Lys31 not found" >> residues.txt
grep " HIS ${CHAIN}  34" ${PDB_ID}.pdb >> residues.txt || echo "His34 not found" >> residues.txt
grep " TYR ${CHAIN}  83" ${PDB_ID}.pdb >> residues.txt || echo "Tyr83 not found" >> residues.txt
grep " LYS ${CHAIN} 353" ${PDB_ID}.pdb >> residues.txt || echo "Lys353 not found" >> residues.txt

# Activate conda environment
conda activate docking

# Clean PDB to isolate ACE2 (chain A), remove solvent, add hydrogens
pymol -c ${PDB_ID}.pdb -d "select chain ${CHAIN}; save ${PDB_ID}_ace2.pdb, sele; remove solvent; h_add; save ${PDB_ID}_clean.pdb, sele"
if [ ! -s ${PDB_ID}_clean.pdb ]; then
    echo "Error: ${PDB_ID}_clean.pdb is empty. Check input PDB and chain."
    exit 1
fi

# Convert cleaned PDB to PDBQT for receptor
obabel ${PDB_ID}_clean.pdb -O ${PDB_ID}_clean.pdbqt -xr -xh -p 7.4 --verbose 0
ls -l ${PDB_ID}_clean.pdbqt

# Verify active site residues in cleaned receptor
pymol -c ${PDB_ID}_clean.pdb -d "select active_site, resi ${ACTIVE_SITE_RESIDUES} and chain ${CHAIN}; print(len(active_site))" > active_site.txt
cat active_site.txt

# Calculate active site center
pymol -c ${PDB_ID}_clean.pdb -d "select active_site, resi ${ACTIVE_SITE_RESIDUES} and chain ${CHAIN}; run https://raw.githubusercontent.com/Pymol-Scripts/Pymol-scripts/master/com.py; get_com active_site" > box_center.txt
cat box_center.txt
# Set box center (adjust based on box_center.txt)
set CENTER_X=-32.8
set CENTER_Y=28.7
set CENTER_Z=17.5

# Generate ligand PDB from SMILES using RDKit
echo "from rdkit import Chem
from rdkit.Chem import AllChem
mol = Chem.MolFromSmiles('${SMILES}')
if mol:
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=88)
    AllChem.UFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, 'remdesivir_rdkit.pdb')
    print('Ligand PDB generated')" > generate_remdesivir_pdb.py
python generate_remdesivir_pdb.py

# Align ligand to active site center
echo "from rdkit import Chem
mol = Chem.MolFromPDBFile('remdesivir_rdkit.pdb')
if mol:
    conf = mol.GetConformer()
    x_sum, y_sum, z_sum, n = 0, 0, 0, 0
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        x_sum += pos.x
        y_sum += pos.y
        z_sum += pos.z
        n += 1
    x_center = x_sum / n
    y_center = y_sum / n
    z_center = z_sum / n
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (pos.x - x_center + ${CENTER_X}, pos.y - y_center + ${CENTER_Y}, pos.z - z_center + ${CENTER_Z}))
    Chem.MolToPDBFile(mol, 'remdesivir_aligned.pdb')
    print('Ligand aligned to (${CENTER_X}, ${CENTER_Y}, ${CENTER_Z})')" > align_remdesivir.py
python align_remdesivir.py

# Convert aligned PDB to PDBQT and remove REMARK tags
obabel remdesivir_aligned.pdb -O remdesivir_aligned.pdbqt -xl -xh -p 7.4 --verbose 0
sed -i '/^REMARK/d' remdesivir_aligned.pdbqt
cat -n remdesivir_aligned.pdbqt | head -n 20

# Deactivate conda environment
conda deactivate

# Perform docking with GNINA
docker run --gpus all -v $(pwd):/work -w /work gnina/gnina:latest gnina \
  --receptor ${PDB_ID}_clean.pdbqt \
  --ligand remdesivir_aligned.pdbqt \
  --center_x ${CENTER_X} --center_y ${CENTER_Y} --center_z ${CENTER_Z} \
  --size_x 30 --size_y 30 --size_z 30 \
  --out test_gnina_specific.pdbqt \
  --cnn_scoring rescore \
  --seed 88 \
  --exhaustiveness 32

# Verify docking output
if [ -f "test_gnina_specific.pdbqt" ]; then
    echo "Success: test_gnina_specific.pdbqt generated."
    ls -la test_gnina_specific.pdbqt
else
    echo "Error: Docking output not generated. Check logs."
    exit 1
fi

# Activate conda environment for visualization
conda activate docking

# Visualize docking results in PyMOL
pymol -c test_gnina_specific.pdbqt ${PDB_ID}_clean.pdbqt -d "select active_site, resi ${ACTIVE_SITE_RESIDUES} and chain ${CHAIN}; show sticks, active_site; label active_site, resi + resn; zoom active_site"

# Deactivate conda environment
conda deactivate
