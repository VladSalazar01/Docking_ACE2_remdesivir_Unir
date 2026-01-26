#!/usr/bin/env python3
"""
EJECUCIÓN DE DOCKING BATCH CON GNINA (DOCKER)
Carpeta de trabajo: G:\Cribado_VC_ACE_2\validacion_ACE1
"""

import subprocess
import sys
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(r"G:\Cribado_VC_ACE_2\validacion_ACE1")
RECEPTOR = WORK_DIR / "receptor.pdb"
LIGANDS = WORK_DIR / "subset" / "library_combined.sdf"
OUTPUT_DIR = WORK_DIR / "results"
OUTPUT = OUTPUT_DIR / "docked_results.sdf"

CENTER = (45.2, 45.6, 45.4)
SIZE = (26, 26, 26)
EXHAUSTIVENESS = 64

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

def verificar_requisitos():
    log("Verificando requisitos...")
    if not RECEPTOR.exists():
        log(f"ERROR: Receptor no encontrado: {RECEPTOR}")
        sys.exit(1)
    log(f"  Receptor: {RECEPTOR}")
    if not LIGANDS.exists():
        log(f"ERROR: Biblioteca no encontrada: {LIGANDS}")
        log("  Ejecutar primero: python build_validation_subset.py")
        sys.exit(1)
    log(f"  Ligandos: {LIGANDS}")
    with open(LIGANDS, 'r') as f:
        n_ligands = f.read().count('$$$$')
    log(f"  Total ligandos: {n_ligands}")
    tiempo_hrs = n_ligands * 25 / 3600
    log(f"  Tiempo estimado: {tiempo_hrs:.1f} horas")
    return n_ligands

def ejecutar_gnina():
    OUTPUT_DIR.mkdir(exist_ok=True)
    work_dir_docker = str(WORK_DIR).replace('\\', '/')
    cmd = [
        "docker", "run", "--rm", "--gpus", "all",
        "-v", f"{work_dir_docker}:/data",
        "gnina/gnina",
        "gnina",
        "-r", "/data/receptor.pdb",
        "-l", "/data/subset/library_combined.sdf",
        "-o", "/data/results/docked_results.sdf",
        "--center_x", str(CENTER[0]),
        "--center_y", str(CENTER[1]),
        "--center_z", str(CENTER[2]),
        "--size_x", str(SIZE[0]),
        "--size_y", str(SIZE[1]),
        "--size_z", str(SIZE[2]),
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--cnn_scoring", "rescore",
        "--num_modes", "1",
        "--device", "0"
    ]
    log("")
    log("Comando Docker:")
    print("  " + " ".join(cmd[:7]))
    print("  " + " ".join(cmd[7:]))
    log("")
    log("Iniciando docking...")
    log("="*60)
    result = subprocess.run(cmd)
    log("="*60)
    if result.returncode == 0:
        log(f"Docking completado: {OUTPUT}")
    else:
        log(f"ERROR: Gnina retorno codigo {result.returncode}")
        sys.exit(1)

def main():
    log("="*60)
    log("VALIDACION RETROSPECTIVA CON GNINA")
    log("TFM: Cribado Virtual ACE2")
    log("="*60)
    n_ligands = verificar_requisitos()
    respuesta = input(f"\nProceder con docking de {n_ligands} ligandos? [s/N]: ")
    if respuesta.lower() != 's':
        log("Cancelado por usuario")
        sys.exit(0)
    ejecutar_gnina()
    log("")
    log("Siguiente paso: python calculate_metrics.py")

if __name__ == "__main__":
    main()