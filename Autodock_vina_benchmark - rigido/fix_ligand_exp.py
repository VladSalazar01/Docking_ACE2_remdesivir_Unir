from rdkit import Chem
from rdkit.Chem import AllChem

print("Leyendo ligand_exp.pdb original...")
mol = Chem.MolFromPDBFile('ligand_exp.pdb', removeHs=False, sanitize=False)

if mol is None:
    print("Error: No se pudo leer ligand_exp.pdb")
    exit(1)

print(f"Átomos originales: {mol.GetNumAtoms()}")

# Contar hidrógenos reales
h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
print(f"Hidrógenos explícitos: {h_count}")
print(f"Átomos pesados: {mol.GetNumAtoms() - h_count}")

# Ver qué hay en la molécula
print("\nPrimeros 10 átomos:")
for i in range(min(10, mol.GetNumAtoms())):
    atom = mol.GetAtomWithIdx(i)
    print(f"  {i}: {atom.GetSymbol()} (num atómico: {atom.GetAtomicNum()})")

# Obtener fragmentos
frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
print(f"\nFragmentos encontrados: {len(frags)}")

for i, frag in enumerate(frags):
    print(f"  Fragmento {i+1}: {frag.GetNumAtoms()} átomos")

if len(frags) > 1:
    # Seleccionar el fragmento más grande
    mol_clean = max(frags, key=lambda m: m.GetNumAtoms())
    print(f"\n✓ Seleccionado fragmento más grande: {mol_clean.GetNumAtoms()} átomos")
else:
    mol_clean = mol

# Intentar sanitizar
try:
    Chem.SanitizeMol(mol_clean)
    print("✓ Molécula sanitizada")
except Exception as e:
    print(f"⚠ No se pudo sanitizar completamente: {e}")

# Remover cualquier hidrógeno
mol_final = Chem.RemoveHs(mol_clean, implicitOnly=False, updateExplicitCount=True)

print(f"\nÁtomos finales: {mol_final.GetNumAtoms()}")

# Guardar
Chem.MolToPDBFile(mol_final, 'ligand_exp_clean.pdb')
print("✓ Guardado como ligand_exp_clean.pdb")

# Verificar
test = Chem.MolFromPDBFile('ligand_exp_clean.pdb', removeHs=False)
if test:
    print(f"✓ Verificación: {test.GetNumAtoms()} átomos")
    
    # Comparar con output_1.pdb
    probe = Chem.MolFromPDBFile('output_1.pdb', removeHs=False)
    if probe:
        print(f"\nComparación:")
        print(f"  ligand_exp_clean.pdb: {test.GetNumAtoms()} átomos")
        print(f"  output_1.pdb:         {probe.GetNumAtoms()} átomos")
        
        if test.GetNumAtoms() == probe.GetNumAtoms():
            print("\n✓✓✓ AMBAS MOLÉCULAS TIENEN EL MISMO NÚMERO DE ÁTOMOS ✓✓✓")
        else:
            print("\n⚠ Todavía hay diferencia en el número de átomos")
else:
    print("⚠ No se pudo verificar el archivo")