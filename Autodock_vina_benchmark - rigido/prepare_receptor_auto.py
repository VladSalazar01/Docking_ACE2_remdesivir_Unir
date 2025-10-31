import subprocess
import os

OBABEL = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"

print("Paso 1: Extrayendo proteína desde 9FMM.cif...")
subprocess.run([
    OBABEL,
    "-i", "cif", "9FMM.cif",
    "-o", "pdb",
    "-O", "9FMM_protein.pdb",
    "--filter", "resname!=A1IDX and resname!=NAG and resname!=EDO and resname!=CL and resname!=NA and resname!=HOH"
], check=True)

if not os.path.exists('9FMM_protein.pdb'):
    print("Error: No se pudo crear 9FMM_protein.pdb")
    exit(1)

print("Paso 2: Convirtiendo a PDBQT con cargas Gasteiger...")
subprocess.run([
    OBABEL,
    "-i", "pdb", "9FMM_protein.pdb",
    "-o", "pdbqt",
    "-O", "receptor.pdbqt",
    "-xr",  # Receptor rígido
    "--partialcharge", "gasteiger"
], check=True)

if os.path.exists('receptor.pdbqt'):
    # Contar líneas
    with open('receptor.pdbqt', 'r') as f:
        lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')]
        print(f"\n✓ Receptor preparado: receptor.pdbqt")
        print(f"  Total de átomos: {len(atom_lines)}")
else:
    print("Error: No se pudo crear receptor.pdbqt")