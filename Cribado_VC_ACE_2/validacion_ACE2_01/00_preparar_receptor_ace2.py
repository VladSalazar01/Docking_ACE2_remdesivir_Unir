#!/usr/bin/env python3
"""
00_preparar_receptor_ace2.py
============================
Preparación del receptor ACE2 (PDB 9FMM) para docking con Gnina.

PROTOCOLO DE PREPARACIÓN:
1. Descargar estructura PDB 9FMM
2. Eliminar aguas, iones y ligando cristalográfico
3. Seleccionar cadena A (monómero funcional)
4. Agregar hidrógenos polares
5. Convertir a formato PDBQT

JUSTIFICACIÓN PARÁMETROS:
- Sin aguas: Gnina no parametriza aguas explícitas
- Solo cadena A: ACE2 es monómero funcional
- Hidrógenos polares: Necesarios para cálculo de puentes de hidrógeno
- Zinc catalítico: MANTENER (esencial para actividad enzimática)

CENTRO DE DOCKING (basado en ligando F-MLN-4760):
- Center: (42.0, 7.0, 23.0) Å
- Size: 26×26×26 Å (cubre sitio activo completo)

Requisitos:
- Open Babel (obabel) instalado
- Opcionalmente: PyMOL o UCSF Chimera para inspección

Autor: Script generado para TFM - Validación Gnina ACE2
Fecha: 2026-01-26
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


# ============================================================================
# CONFIGURACIÓN
# ============================================================================

CONFIG = {
    "pdb_id": "9FMM",
    "cadena": "A",
    "output_name": "9FMM_receptor.pdbqt",
    
    # Centro del sitio activo (coordenadas del ligando F-MLN-4760)
    "center_x": 42.0,
    "center_y": 7.0,
    "center_z": 23.0,
    
    # Tamaño de la caja de docking
    "size_x": 26,
    "size_y": 26,
    "size_z": 26,
}


# ============================================================================
# FUNCIONES
# ============================================================================

def verificar_dependencias():
    """Verifica que Open Babel esté instalado."""
    try:
        result = subprocess.run(["obabel", "-V"], capture_output=True, text=True)
        print(f"✓ Open Babel detectado: {result.stdout.strip()}")
        return True
    except FileNotFoundError:
        print("✗ Open Babel no encontrado")
        print("\nInstalación:")
        print("  Windows: conda install -c conda-forge openbabel")
        print("  Linux:   sudo apt install openbabel")
        return False


def descargar_pdb(pdb_id: str, output_path: str) -> bool:
    """Busca estructura local o descarga desde RCSB."""
    import urllib.request
    
    # Primero buscar archivo local (puede ser .cif o .pdb)
    local_cif = Path(f"{pdb_id}.cif")
    local_pdb = Path(f"{pdb_id}.pdb")
    
    if local_cif.exists():
        print(f"✓ Archivo local encontrado: {local_cif}")
        # Convertir CIF a PDB usando Open Babel
        cmd = ["obabel", str(local_cif), "-O", output_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ Convertido CIF → PDB: {output_path}")
            return True
        else:
            print(f"⚠ Error convirtiendo CIF: {result.stderr}")
            # Intentar procesamiento manual
            return convertir_cif_a_pdb_manual(str(local_cif), output_path)
    
    if local_pdb.exists():
        print(f"✓ Archivo local encontrado: {local_pdb}")
        # Copiar al path de salida
        shutil.copy(local_pdb, output_path)
        return True
    
    # Si no hay archivo local, intentar descarga
    # Primero intentar .pdb, luego .cif
    for ext, fmt in [(".pdb", "pdb"), (".cif", "cif")]:
        url = f"https://files.rcsb.org/download/{pdb_id}{ext}"
        try:
            print(f"Descargando {pdb_id}{ext} desde RCSB...")
            urllib.request.urlretrieve(url, f"{pdb_id}{ext}")
            
            if ext == ".cif":
                # Convertir a PDB
                cmd = ["obabel", f"{pdb_id}{ext}", "-O", output_path]
                subprocess.run(cmd, capture_output=True)
            else:
                shutil.copy(f"{pdb_id}{ext}", output_path)
            
            print(f"✓ Descargado: {output_path}")
            return True
        except Exception as e:
            print(f"  {ext}: {e}")
            continue
    
    print(f"✗ No se pudo obtener {pdb_id}")
    return False


def convertir_cif_a_pdb_manual(cif_path: str, pdb_path: str) -> bool:
    """
    Conversión manual básica de mmCIF a PDB.
    Extrae coordenadas ATOM/HETATM del archivo CIF.
    """
    print("  Intentando conversión manual CIF → PDB...")
    
    lineas_pdb = []
    
    try:
        with open(cif_path, "r") as f:
            en_atom_site = False
            columnas = {}
            
            for linea in f:
                linea = linea.strip()
                
                # Detectar sección _atom_site
                if linea.startswith("_atom_site."):
                    en_atom_site = True
                    col_name = linea.split(".")[1].split()[0]
                    columnas[col_name] = len(columnas)
                    continue
                
                if en_atom_site and linea.startswith("_"):
                    en_atom_site = False
                    continue
                
                if en_atom_site and linea and not linea.startswith("#"):
                    # Parsear línea de átomos
                    partes = linea.split()
                    if len(partes) < 10:
                        continue
                    
                    try:
                        group = partes[columnas.get("group_PDB", 0)]
                        atom_id = partes[columnas.get("id", 1)]
                        atom_name = partes[columnas.get("label_atom_id", 3)]
                        res_name = partes[columnas.get("label_comp_id", 5)]
                        chain = partes[columnas.get("label_asym_id", 6)]
                        res_seq = partes[columnas.get("label_seq_id", 8)]
                        x = float(partes[columnas.get("Cartn_x", 10)])
                        y = float(partes[columnas.get("Cartn_y", 11)])
                        z = float(partes[columnas.get("Cartn_z", 12)])
                        element = partes[columnas.get("type_symbol", 2)]
                        
                        # Formatear línea PDB
                        record = "ATOM  " if group == "ATOM" else "HETATM"
                        pdb_line = f"{record}{int(atom_id):5d} {atom_name:<4s} {res_name:>3s} {chain:1s}{int(res_seq) if res_seq != '.' else 0:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}\n"
                        lineas_pdb.append(pdb_line)
                    except (ValueError, IndexError, KeyError):
                        continue
        
        if lineas_pdb:
            lineas_pdb.append("END\n")
            with open(pdb_path, "w") as f:
                f.writelines(lineas_pdb)
            print(f"  ✓ Conversión manual exitosa: {len(lineas_pdb)-1} átomos")
            return True
        
    except Exception as e:
        print(f"  ✗ Error en conversión manual: {e}")
    
    return False


def limpiar_pdb(input_pdb: str, output_pdb: str, cadena: str = "A"):
    """
    Limpia el PDB: elimina aguas, heteroátomos no esenciales, 
    selecciona cadena y mantiene zinc catalítico.
    """
    
    lineas_limpias = []
    zinc_encontrado = False
    
    with open(input_pdb, "r") as f:
        for linea in f:
            # Saltar aguas
            if linea.startswith(("HETATM", "ATOM")):
                res_name = linea[17:20].strip()
                chain = linea[21]
                
                # Saltar aguas
                if res_name == "HOH":
                    continue
                
                # Saltar ligando cristalográfico (F-MLN-4760 = 9FM)
                if res_name == "9FM":
                    continue
                
                # Solo cadena especificada
                if chain != cadena:
                    continue
                
                # Mantener zinc catalítico (esencial para ACE2)
                if res_name == "ZN":
                    zinc_encontrado = True
                    lineas_limpias.append(linea)
                    continue
                
                # Mantener cloruro si está coordinado
                if res_name == "CL":
                    lineas_limpias.append(linea)
                    continue
                
                # Saltar otros heteroátomos (buffer, etc)
                if linea.startswith("HETATM") and res_name not in ["ZN", "CL"]:
                    continue
                
                lineas_limpias.append(linea)
            
            elif linea.startswith("END"):
                lineas_limpias.append(linea)
                break
    
    with open(output_pdb, "w") as f:
        f.writelines(lineas_limpias)
    
    n_atoms = sum(1 for l in lineas_limpias if l.startswith("ATOM"))
    print(f"✓ PDB limpio: {n_atoms} átomos, cadena {cadena}")
    
    if zinc_encontrado:
        print("✓ Zinc catalítico preservado (esencial para ACE2)")
    else:
        print("⚠ Zinc no encontrado - verificar estructura")
    
    return output_pdb


def agregar_hidrogenos(input_pdb: str, output_pdb: str):
    """Agrega hidrógenos polares usando Open Babel."""
    
    cmd = [
        "obabel", input_pdb,
        "-O", output_pdb,
        "-h",  # Agregar hidrógenos
        "-p", "7.4",  # pH fisiológico
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"✓ Hidrógenos agregados (pH 7.4)")
        return True
    else:
        print(f"✗ Error agregando hidrógenos: {result.stderr}")
        return False


def convertir_pdbqt(input_pdb: str, output_pdbqt: str):
    """Convierte PDB a PDBQT usando Open Babel."""
    
    cmd = [
        "obabel", input_pdb,
        "-O", output_pdbqt,
        "-xr",  # Receptor (rígido)
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0 and os.path.exists(output_pdbqt):
        size = os.path.getsize(output_pdbqt)
        print(f"✓ PDBQT generado: {output_pdbqt} ({size:,} bytes)")
        return True
    else:
        print(f"✗ Error convirtiendo a PDBQT: {result.stderr}")
        return False


def generar_archivo_config(output_dir: str):
    """Genera archivo de configuración para Gnina/Vina."""
    
    config_path = Path(output_dir) / "config_docking.txt"
    
    contenido = f"""# Configuración de docking ACE2 (PDB 9FMM)
# Generado automáticamente por 00_preparar_receptor_ace2.py

# Receptor
receptor = {CONFIG['output_name']}

# Centro del sitio activo (basado en F-MLN-4760)
center_x = {CONFIG['center_x']}
center_y = {CONFIG['center_y']}
center_z = {CONFIG['center_z']}

# Tamaño de la caja de búsqueda
size_x = {CONFIG['size_x']}
size_y = {CONFIG['size_y']}
size_z = {CONFIG['size_z']}

# Parámetros de búsqueda
exhaustiveness = 64
num_modes = 1

# Notas:
# - El sitio activo contiene Zn(II) catalitico
# - Centro basado en coordenadas del ligando cristalográfico
# - Caja de 26Å cubre completamente el bolsillo de unión
"""
    
    with open(config_path, "w", encoding="utf-8") as f:
        f.write(contenido)
    
    print(f"✓ Configuración guardada: {config_path}")
    return config_path


def main():
    """Flujo principal de preparación del receptor."""
    
    print("=" * 70)
    print("PREPARACIÓN DEL RECEPTOR ACE2 (PDB 9FMM)")
    print("=" * 70)
    
    # Verificar dependencias
    if not verificar_dependencias():
        sys.exit(1)
    
    # Crear directorio de trabajo
    work_dir = Path(".")
    
    # Archivos intermedios
    pdb_original = work_dir / f"{CONFIG['pdb_id']}.pdb"
    pdb_limpio = work_dir / f"{CONFIG['pdb_id']}_limpio.pdb"
    pdb_hidrogenado = work_dir / f"{CONFIG['pdb_id']}_H.pdb"
    pdbqt_final = work_dir / CONFIG['output_name']
    
    # También verificar si existe .cif
    cif_local = work_dir / f"{CONFIG['pdb_id']}.cif"
    
    print("-" * 50)
    print("PASO 1: Obtener estructura PDB")
    print("-" * 50)
    
    if pdb_original.exists():
        print(f"✓ PDB ya existe: {pdb_original}")
    elif cif_local.exists():
        print(f"✓ CIF local encontrado: {cif_local}")
        if not descargar_pdb(CONFIG['pdb_id'], str(pdb_original)):
            sys.exit(1)
    else:
        if not descargar_pdb(CONFIG['pdb_id'], str(pdb_original)):
            sys.exit(1)
    
    print("-" * 50)
    print("PASO 2: Limpiar estructura")
    print("-" * 50)
    print("  - Eliminando aguas (HOH)")
    print("  - Eliminando ligando cristalográfico (9FM)")
    print("  - Seleccionando cadena A")
    print("  - Preservando Zn catalítico")
    
    limpiar_pdb(str(pdb_original), str(pdb_limpio), CONFIG['cadena'])
    
    print("-" * 50)
    print("PASO 3: Agregar hidrógenos polares")
    print("-" * 50)
    
    if not agregar_hidrogenos(str(pdb_limpio), str(pdb_hidrogenado)):
        # Alternativa: usar PDB limpio sin hidrógenos extras
        print("⚠ Usando PDB sin hidrógenos adicionales")
        pdb_hidrogenado = pdb_limpio
    
    print("-" * 50)
    print("PASO 4: Convertir a PDBQT")
    print("-" * 50)
    
    if not convertir_pdbqt(str(pdb_hidrogenado), str(pdbqt_final)):
        sys.exit(1)
    
    print("-" * 50)
    print("PASO 5: Generar archivo de configuración")
    print("-" * 50)
    
    generar_archivo_config(str(work_dir))
    
    # Resumen
    print("=" * 70)
    print("PREPARACIÓN COMPLETADA")
    print("=" * 70)
    print(f"""
Archivos generados:
  - {pdbqt_final} (receptor para docking)
  - config_docking.txt (parámetros de caja)

Parámetros del sitio activo:
  Centro: ({CONFIG['center_x']}, {CONFIG['center_y']}, {CONFIG['center_z']}) Å
  Tamaño: {CONFIG['size_x']}×{CONFIG['size_y']}×{CONFIG['size_z']} Å

VERIFICACIÓN RECOMENDADA:
  1. Abrir {pdbqt_final} en PyMOL/Chimera
  2. Verificar que el Zn catalítico esté presente
  3. Verificar que no haya aguas ni ligando
  4. Comprobar que la caja cubra el sitio activo

Comando de verificación en PyMOL:
  load {pdbqt_final}
  show cartoon
  show spheres, resn ZN
  pseudoatom center, pos=[{CONFIG['center_x']},{CONFIG['center_y']},{CONFIG['center_z']}]
  show spheres, center
""")
    
    # Limpiar archivos intermedios (opcional)
    # pdb_limpio.unlink()
    # pdb_hidrogenado.unlink()
    
    return str(pdbqt_final)


if __name__ == "__main__":
    receptor = main()