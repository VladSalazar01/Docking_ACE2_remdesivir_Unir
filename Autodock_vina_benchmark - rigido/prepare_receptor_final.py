import subprocess
import os

OBABEL = r"C:\Program Files\OpenBabel-3.1.1\obabel.exe"

print("Convirtiendo 9FMM_protein.pdb a receptor.pdbqt...")

if not os.path.exists('9FMM_protein.pdb'):
    print("Error: 9FMM_protein.pdb no encontrado")
    exit(1)

# Convertir a PDBQT
result = subprocess.run([
    OBABEL,
    "-i", "pdb", "9FMM_protein.pdb",
    "-o", "pdbqt",
    "-O", "receptor.pdbqt",
    "-xr",  # Receptor rígido
    "--partialcharge", "gasteiger"
], capture_output=True, text=True)

print(result.stdout)
if result.stderr and "converted" not in result.stderr:
    print("Warnings:", result.stderr)

# Verificar resultado
if os.path.exists('receptor.pdbqt'):
    with open('receptor.pdbqt', 'r') as f:
        lines = f.readlines()
        atom_lines = [l for l in lines if l.startswith('ATOM')]
    
    if len(atom_lines) > 0:
        print(f"\n✓ Receptor preparado exitosamente")
        print(f"  Archivo: receptor.pdbqt")
        print(f"  Átomos: {len(atom_lines)}")
    else:
        print("\n✗ Error: receptor.pdbqt está vacío")
        print("Verifica que 9FMM_protein.pdb contenga átomos ATOM")
else:
    print("\n✗ Error: No se pudo crear receptor.pdbqt")