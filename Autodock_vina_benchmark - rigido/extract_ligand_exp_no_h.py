# extract_ligand_exp_no_h.py
from pymol import cmd

cmd.load('9FMM.cif')
cmd.select('ligand', 'resn A1IDX')

# Remover hidrógenos ANTES de guardar
cmd.remove('ligand and elem H')

# Guardar solo átomos pesados
cmd.save('ligand_exp.pdb', 'ligand')

atoms = cmd.count_atoms('ligand')
print(f"✓ Pose experimental extraída SIN hidrógenos: ligand_exp.pdb ({atoms} átomos)")

cmd.quit()