from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem

# Leer ligando
mol = Chem.MolFromPDBFile('ligand_A1IDX.pdb', removeHs=False)

if mol is None:
    print("Error: No se pudo leer ligand_A1IDX.pdb")
else:
    # Agregar hidrógenos explícitos
    mol = Chem.AddHs(mol, addCoords=True)
    
    # Obtener el fragmento más grande (eliminar aguas/iones)
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:
        print(f"⚠ Encontrados {len(frags)} fragmentos. Seleccionando el más grande...")
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Preparar para docking
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    
    # Escribir PDBQT
    for setup in mol_setups:
        pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        
        if is_ok:
            with open('ligand.pdbqt', 'w') as f:
                f.write(pdbqt_string)
            print("✓ Ligando preparado: ligand.pdbqt")
            print(f"  Átomos: {mol.GetNumAtoms()}")
        else:
            print(f"Error al escribir PDBQT: {error_msg}")