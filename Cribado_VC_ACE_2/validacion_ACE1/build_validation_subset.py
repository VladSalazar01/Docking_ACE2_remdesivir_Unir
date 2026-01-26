#!/usr/bin/env python3
"""
CONSTRUCCIÓN DE SUBCONJUNTO DE VALIDACIÓN RETROSPECTIVA
Carpeta de trabajo: G:\Cribado_VC_ACE_2\validacion_ACE1
"""

import os
import sys
import gzip
import random
import subprocess
import shutil
from datetime import datetime
from pathlib import Path

WORK_DIR = Path(r"G:\Cribado_VC_ACE_2\validacion_ACE1")
SEED = 42
N_ACTIVES = 50
N_DECOYS = 1000

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

def parse_mol2(filename):
    molecules = []
    current_mol = []
    if str(filename).endswith('.gz'):
        f = gzip.open(filename, 'rt', encoding='utf-8', errors='ignore')
    else:
        f = open(filename, 'r', encoding='utf-8', errors='ignore')
    try:
        for line in f:
            if line.startswith('@<TRIPOS>MOLECULE'):
                if current_mol:
                    molecules.append(''.join(current_mol))
                current_mol = [line]
            elif current_mol:
                current_mol.append(line)
        if current_mol:
            molecules.append(''.join(current_mol))
    finally:
        f.close()
    return molecules

def extract_names_from_sdf(filename):
    names = []
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        current_name = None
        for line in f:
            if current_name is None and line.strip() and not line.startswith(' '):
                current_name = line.strip()
            if line.startswith('$$$$'):
                if current_name:
                    names.append(current_name)
                current_name = None
    return names

def paso1_verificar():
    log("="*60)
    log("PASO 1: VERIFICACION DE ARCHIVOS")
    log("="*60)
    os.chdir(WORK_DIR)
    log(f"Directorio de trabajo: {WORK_DIR}")
    actives_file = None
    decoys_file = None
    for pattern in ['actives_final.mol2.gz', 'actives_final.mol2', 'actives_final_mol2.gz']:
        if Path(pattern).exists():
            actives_file = pattern
            break
    for pattern in ['decoys_final.mol2.gz', 'decoys_final.mol2', 'decoys_final_mol2.gz']:
        if Path(pattern).exists():
            decoys_file = pattern
            break
    tar_file = Path('dude_ace1_validation_package.tar.gz')
    if tar_file.exists() and (actives_file is None or decoys_file is None):
        log("Extrayendo archivos de dude_ace1_validation_package.tar.gz...")
        import tarfile
        with tarfile.open(tar_file, 'r:gz') as tar:
            tar.extractall()
        src_dir = Path('dude_ace1_validation')
        if src_dir.exists():
            for f in src_dir.glob('*.mol2*'):
                shutil.copy(f, WORK_DIR)
            actives_file = 'actives_final.mol2.gz' if Path('actives_final.mol2.gz').exists() else 'actives_final.mol2'
            decoys_file = 'decoys_final.mol2.gz' if Path('decoys_final.mol2.gz').exists() else 'decoys_final.mol2'
    if actives_file is None:
        log("ERROR: No se encontro archivo de activos")
        sys.exit(1)
    if decoys_file is None:
        log("ERROR: No se encontro archivo de decoys")
        sys.exit(1)
    log(f"  Activos: {actives_file}")
    log(f"  Decoys: {decoys_file}")
    try:
        result = subprocess.run(['obabel', '-V'], capture_output=True, text=True)
        version = result.stdout.split()[2] if result.stdout else 'disponible'
        log(f"  Open Babel: {version}")
    except FileNotFoundError:
        log("ERROR: Open Babel no encontrado")
        sys.exit(1)
    return actives_file, decoys_file

def paso2_subconjunto(actives_file, decoys_file):
    log("")
    log("="*60)
    log("PASO 2: SELECCION DE SUBCONJUNTO ESTRATIFICADO")
    log("="*60)
    random.seed(SEED)
    log(f"Semilla aleatoria: {SEED}")
    log(f"Cargando activos desde {actives_file}...")
    actives = parse_mol2(actives_file)
    log(f"  {len(actives)} activos cargados")
    log(f"Cargando decoys desde {decoys_file}...")
    decoys = parse_mol2(decoys_file)
    log(f"  {len(decoys)} decoys cargados")
    log(f"Seleccionando: {N_ACTIVES} activos + {N_DECOYS} decoys")
    subset_actives = random.sample(actives, min(N_ACTIVES, len(actives)))
    subset_decoys = random.sample(decoys, min(N_DECOYS, len(decoys)))
    ratio = len(subset_decoys) / len(subset_actives)
    log(f"  Ratio decoys:activos = {ratio:.0f}:1")
    return subset_actives, subset_decoys

def paso3_convertir(subset_actives, subset_decoys):
    log("")
    log("="*60)
    log("PASO 3: CONVERSION MOL2 -> SDF")
    log("="*60)
    subset_dir = WORK_DIR / 'subset'
    subset_dir.mkdir(exist_ok=True)
    log(f"Directorio: {subset_dir}")
    actives_mol2 = subset_dir / 'actives_subset.mol2'
    with open(actives_mol2, 'w') as f:
        f.writelines(subset_actives)
    log(f"  actives_subset.mol2: {len(subset_actives)} moleculas")
    decoys_mol2 = subset_dir / 'decoys_subset.mol2'
    with open(decoys_mol2, 'w') as f:
        f.writelines(subset_decoys)
    log(f"  decoys_subset.mol2: {len(subset_decoys)} moleculas")
    log("Convirtiendo con Open Babel...")
    actives_sdf = subset_dir / 'actives_subset.sdf'
    subprocess.run(['obabel', str(actives_mol2), '-O', str(actives_sdf)], capture_output=True)
    log(f"  actives_subset.sdf")
    decoys_sdf = subset_dir / 'decoys_subset.sdf'
    subprocess.run(['obabel', str(decoys_mol2), '-O', str(decoys_sdf)], capture_output=True)
    log(f"  decoys_subset.sdf")
    library_sdf = subset_dir / 'library_combined.sdf'
    with open(library_sdf, 'w') as out:
        with open(actives_sdf, 'r') as f:
            out.write(f.read())
        with open(decoys_sdf, 'r') as f:
            out.write(f.read())
    n_total = len(subset_actives) + len(subset_decoys)
    log(f"  library_combined.sdf: {n_total} moleculas")
    return actives_sdf, decoys_sdf

def paso4_etiquetas(actives_sdf, decoys_sdf):
    log("")
    log("="*60)
    log("PASO 4: GENERACION DE ETIQUETAS")
    log("="*60)
    actives_names = extract_names_from_sdf(actives_sdf)
    decoys_names = extract_names_from_sdf(decoys_sdf)
    labels_file = WORK_DIR / 'subset' / 'labels.csv'
    with open(labels_file, 'w') as f:
        f.write("mol_name,is_active\n")
        for name in actives_names:
            f.write(f"{name},1\n")
        for name in decoys_names:
            f.write(f"{name},0\n")
    log(f"  labels.csv: {len(actives_names) + len(decoys_names)} entradas")
    log(f"    Activos (is_active=1): {len(actives_names)}")
    log(f"    Decoys (is_active=0): {len(decoys_names)}")

def main():
    log("="*60)
    log("CONSTRUCCION DE SUBCONJUNTO DE VALIDACION")
    log("TFM: Cribado Virtual ACE2")
    log("="*60)
    actives_file, decoys_file = paso1_verificar()
    subset_actives, subset_decoys = paso2_subconjunto(actives_file, decoys_file)
    actives_sdf, decoys_sdf = paso3_convertir(subset_actives, subset_decoys)
    paso4_etiquetas(actives_sdf, decoys_sdf)
    log("")
    log("="*60)
    log("COMPLETADO")
    log("="*60)
    n_total = N_ACTIVES + N_DECOYS
    tiempo = n_total * 25 / 3600
    log(f"Subconjunto: {N_ACTIVES} activos + {N_DECOYS} decoys")
    log(f"Tiempo estimado Gnina: {tiempo:.1f} horas")
    log("")
    log("Archivos generados en subset/:")
    log("  - library_combined.sdf")
    log("  - labels.csv")

if __name__ == "__main__":
    main()