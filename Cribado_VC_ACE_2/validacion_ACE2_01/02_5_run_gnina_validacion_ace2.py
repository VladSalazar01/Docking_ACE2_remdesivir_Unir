#!/usr/bin/env python3
"""
run_gnina_validacion_ace2.py
============================
Ejecuta validación retrospectiva de Gnina para ACE2 usando el dataset generado.

Configuración ACE2:
- Receptor: PDB 9FMM (ACE2 humana con F-MLN-4760)
- Centro grid: (42.0, 7.0, 23.0) - calculado desde ligando cristalográfico
- Tamaño box: 26x26x26 Å
- Exhaustiveness: 64
- CNN scoring: rescore

Métricas calculadas:
- AUC-ROC
- Enrichment Factor (EF) a 1%, 5%, 10%
- BEDROC (α=20)

Autor: Claude para TFM - Validación Gnina ACE2
Fecha: 2025-01-25
"""

import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

# Verificar e instalar dependencias
try:
    import numpy as np
except ImportError:
    os.system("pip install numpy --break-system-packages -q")
    import numpy as np

try:
    from sklearn.metrics import roc_auc_score, roc_curve
except ImportError:
    os.system("pip install scikit-learn --break-system-packages -q")
    from sklearn.metrics import roc_auc_score, roc_curve

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    os.system("pip install rdkit --break-system-packages -q")
    from rdkit import Chem
    from rdkit.Chem import AllChem


# ============================================================================
# CONFIGURACIÓN ACE2
# ============================================================================
CONFIG = {
    # Receptor
    "receptor_pdb": "9FMM",  # PDB ID
    "receptor_file": "9FMM_receptor.pdbqt",  # Archivo preparado
    
    # Grid box (calculado desde F-MLN-4760 en 9FMM)
    "center_x": 42.0,
    "center_y": 7.0,
    "center_z": 23.0,
    "size_x": 26,
    "size_y": 26,
    "size_z": 26,
    
    # Parámetros de docking Gnina
    "exhaustiveness": 64,
    "cnn_scoring": "rescore",
    "num_modes": 1,  # Solo mejor pose para validación
    "cpu": 12,
    
    # Directorios
    "input_dir": "dataset_ace2/validacion",
    "output_dir": "resultados_validacion_ace2",
    "docker_image": "gnina/gnina",
    
    # Carpeta de trabajo Windows
    "trabajo_windows": "G:/Cribado_VC_ACE_2/validacion_ACE2_01"
}


def preparar_ligandos_pdbqt(input_smi: str, output_dir: str) -> List[str]:
    """Convierte ligandos de SMI a PDBQT usando RDKit + Open Babel."""
    
    output_path = Path(output_dir) / "ligandos_pdbqt"
    output_path.mkdir(parents=True, exist_ok=True)
    
    ligandos = []
    
    with open(input_smi, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            
            smiles = parts[0]
            name = parts[1]
            
            # Generar 3D con RDKit
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"  WARN: No se pudo parsear {name}")
                continue
            
            mol = Chem.AddHs(mol)
            
            # Generar conformación 3D
            try:
                result = AllChem.EmbedMolecule(mol, randomSeed=42)
                if result == -1:
                    print(f"  WARN: No se pudo generar 3D para {name}")
                    continue
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception as e:
                print(f"  WARN: Error en {name}: {e}")
                continue
            
            # Guardar como SDF temporal
            sdf_file = output_path / f"{name}.sdf"
            pdbqt_file = output_path / f"{name}.pdbqt"
            
            writer = Chem.SDWriter(str(sdf_file))
            writer.write(mol)
            writer.close()
            
            # Convertir a PDBQT con Open Babel
            cmd = f"obabel {sdf_file} -O {pdbqt_file} -h"
            result = subprocess.run(cmd, shell=True, capture_output=True)
            
            if pdbqt_file.exists():
                ligandos.append(str(pdbqt_file))
            
            # Limpiar SDF temporal
            if sdf_file.exists():
                sdf_file.unlink()
    
    print(f"Ligandos preparados: {len(ligandos)}")
    return ligandos


def ejecutar_gnina_docker(receptor: str, ligando: str, output_file: str) -> dict:
    """Ejecuta Gnina en Docker y retorna resultados."""
    
    # Construir comando Docker
    work_dir = Path(ligando).parent.parent.absolute()
    
    cmd = [
        "docker", "run", "--rm", "--gpus", "all",
        "-v", f"{work_dir}:/data",
        CONFIG["docker_image"],
        "gnina",
        "--receptor", f"/data/{Path(receptor).name}",
        "--ligand", f"/data/ligandos_pdbqt/{Path(ligando).name}",
        "--out", f"/data/salidas/{Path(output_file).name}",
        "--center_x", str(CONFIG["center_x"]),
        "--center_y", str(CONFIG["center_y"]),
        "--center_z", str(CONFIG["center_z"]),
        "--size_x", str(CONFIG["size_x"]),
        "--size_y", str(CONFIG["size_y"]),
        "--size_z", str(CONFIG["size_z"]),
        "--exhaustiveness", str(CONFIG["exhaustiveness"]),
        "--cnn_scoring", CONFIG["cnn_scoring"],
        "--num_modes", str(CONFIG["num_modes"]),
        "--cpu", str(CONFIG["cpu"])
    ]
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 min timeout por ligando
        )
        
        # Parsear salida para obtener scores
        return parsear_salida_gnina(output_file, result.stdout)
        
    except subprocess.TimeoutExpired:
        return {"error": "timeout"}
    except Exception as e:
        return {"error": str(e)}


def parsear_salida_gnina(output_file: str, stdout: str) -> dict:
    """Parsea la salida de Gnina para extraer scores."""
    
    resultado = {
        "cnn_score": None,
        "cnn_affinity": None,
        "vina_affinity": None
    }
    
    # Intentar leer del archivo de salida PDBQT
    try:
        if Path(output_file).exists():
            with open(output_file, "r") as f:
                for line in f:
                    if line.startswith("REMARK"):
                        parts = line.split()
                        if "CNN_VS" in line or "CNNscore" in line:
                            resultado["cnn_score"] = float(parts[-1])
                        elif "CNNaffinity" in line:
                            resultado["cnn_affinity"] = float(parts[-1])
                        elif "VINA" in line and "affinity" in line.lower():
                            resultado["vina_affinity"] = float(parts[-1])
    except Exception:
        pass
    
    # Intentar parsear stdout si no hay archivo
    if resultado["cnn_score"] is None and stdout:
        for line in stdout.split("\n"):
            if "CNNscore" in line or "CNN_VS" in line:
                try:
                    resultado["cnn_score"] = float(line.split()[-1])
                except ValueError:
                    pass
    
    return resultado


def calcular_bedroc(y_true: np.ndarray, y_scores: np.ndarray, alpha: float = 20.0) -> float:
    """
    Calcula BEDROC (Boltzmann-Enhanced Discrimination of ROC).
    
    BEDROC enfatiza el enriquecimiento temprano, crucial para virtual screening.
    alpha=20 corresponde aproximadamente a top 8% del dataset.
    """
    n = len(y_true)
    n_actives = np.sum(y_true)
    
    if n_actives == 0 or n_actives == n:
        return 0.0
    
    # Ordenar por scores descendentes
    sorted_indices = np.argsort(-y_scores)
    sorted_labels = y_true[sorted_indices]
    
    # Calcular BEDROC
    positions = np.where(sorted_labels == 1)[0] + 1  # 1-indexed
    s = np.sum(np.exp(-alpha * positions / n))
    
    # Normalización
    Ra = n_actives / n
    Ri = 1 - np.exp(-alpha)
    Rmin = (Ra / n) * (1 - np.exp(-alpha)) / (np.exp(alpha/n) - 1)
    Rmax = (1 - np.exp(-alpha * Ra)) / (Ra * (1 - np.exp(-alpha)))
    
    bedroc = (s / n_actives - Rmin) / (Rmax - Rmin)
    
    return max(0.0, min(1.0, bedroc))


def calcular_enrichment_factor(y_true: np.ndarray, y_scores: np.ndarray, 
                                porcentaje: float) -> float:
    """
    Calcula el Factor de Enriquecimiento (EF) a un porcentaje dado.
    
    EF = (Activos en top X%) / (Activos esperados en top X%)
    """
    n = len(y_true)
    n_actives = np.sum(y_true)
    
    if n_actives == 0:
        return 0.0
    
    # Top X% del dataset
    top_n = max(1, int(n * porcentaje / 100))
    
    # Ordenar por scores descendentes
    sorted_indices = np.argsort(-y_scores)
    top_labels = y_true[sorted_indices[:top_n]]
    
    # Activos encontrados en top X%
    activos_encontrados = np.sum(top_labels)
    
    # Activos esperados por azar
    activos_esperados = n_actives * porcentaje / 100
    
    if activos_esperados == 0:
        return 0.0
    
    ef = activos_encontrados / activos_esperados
    
    return ef


def calcular_metricas(resultados: List[dict]) -> dict:
    """Calcula todas las métricas de validación."""
    
    # Extraer labels y scores
    labels = np.array([r["label"] for r in resultados])
    
    # Usar CNN score como métrica principal (mayor es mejor)
    scores = np.array([r.get("cnn_score", 0) or 0 for r in resultados])
    
    # Invertir scores de afinidad si se usan (menor es mejor -> mayor es mejor)
    if np.mean(scores) < 0:
        scores = -scores
    
    metricas = {}
    
    # AUC-ROC
    try:
        metricas["auc_roc"] = roc_auc_score(labels, scores)
    except ValueError:
        metricas["auc_roc"] = 0.5
    
    # Enrichment Factors
    metricas["ef_1"] = calcular_enrichment_factor(labels, scores, 1.0)
    metricas["ef_5"] = calcular_enrichment_factor(labels, scores, 5.0)
    metricas["ef_10"] = calcular_enrichment_factor(labels, scores, 10.0)
    
    # BEDROC
    metricas["bedroc_20"] = calcular_bedroc(labels, scores, alpha=20.0)
    
    # Estadísticas adicionales
    metricas["n_activos"] = int(np.sum(labels))
    metricas["n_decoys"] = int(len(labels) - np.sum(labels))
    metricas["ratio"] = metricas["n_decoys"] / max(metricas["n_activos"], 1)
    
    return metricas


def generar_reporte(metricas: dict, output_dir: str):
    """Genera reporte de validación en formato TFM."""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Archivo de texto
    reporte_file = output_path / "reporte_validacion_ace2.txt"
    
    with open(reporte_file, "w", encoding="utf-8") as f:
        f.write("="*70 + "\n")
        f.write("REPORTE DE VALIDACIÓN RETROSPECTIVA - ACE2\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Receptor: PDB {CONFIG['receptor_pdb']}\n")
        f.write(f"Programa: Gnina (CNN scoring: {CONFIG['cnn_scoring']})\n\n")
        
        f.write("-"*70 + "\n")
        f.write("CONFIGURACIÓN DEL GRID\n")
        f.write("-"*70 + "\n")
        f.write(f"Centro: ({CONFIG['center_x']}, {CONFIG['center_y']}, {CONFIG['center_z']}) Å\n")
        f.write(f"Tamaño: {CONFIG['size_x']}×{CONFIG['size_y']}×{CONFIG['size_z']} Å\n")
        f.write(f"Exhaustiveness: {CONFIG['exhaustiveness']}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("COMPOSICIÓN DEL DATASET\n")
        f.write("-"*70 + "\n")
        f.write(f"Activos: {metricas['n_activos']}\n")
        f.write(f"Decoys:  {metricas['n_decoys']}\n")
        f.write(f"Ratio:   1:{metricas['ratio']:.1f}\n\n")
        
        f.write("-"*70 + "\n")
        f.write("MÉTRICAS DE VALIDACIÓN\n")
        f.write("-"*70 + "\n")
        f.write(f"{'Métrica':<20} {'Valor':<10} {'Umbral aceptable':<20}\n")
        f.write("-"*50 + "\n")
        f.write(f"{'AUC-ROC':<20} {metricas['auc_roc']:.3f}{'':>4} {'≥0.70':<20}\n")
        f.write(f"{'EF 1%':<20} {metricas['ef_1']:.2f}{'':>5} {'≥10':<20}\n")
        f.write(f"{'EF 5%':<20} {metricas['ef_5']:.2f}{'':>5} {'≥5':<20}\n")
        f.write(f"{'EF 10%':<20} {metricas['ef_10']:.2f}{'':>5} {'≥3':<20}\n")
        f.write(f"{'BEDROC (α=20)':<20} {metricas['bedroc_20']:.3f}{'':>4} {'≥0.30':<20}\n\n")
        
        # Evaluación
        f.write("-"*70 + "\n")
        f.write("EVALUACIÓN\n")
        f.write("-"*70 + "\n")
        
        criterios_cumplidos = 0
        total_criterios = 5
        
        if metricas["auc_roc"] >= 0.70:
            criterios_cumplidos += 1
            f.write("✓ AUC-ROC ≥ 0.70: CUMPLE\n")
        else:
            f.write(f"✗ AUC-ROC ≥ 0.70: NO CUMPLE ({metricas['auc_roc']:.3f})\n")
        
        if metricas["ef_1"] >= 10:
            criterios_cumplidos += 1
            f.write("✓ EF 1% ≥ 10: CUMPLE\n")
        else:
            f.write(f"✗ EF 1% ≥ 10: NO CUMPLE ({metricas['ef_1']:.2f})\n")
        
        if metricas["ef_5"] >= 5:
            criterios_cumplidos += 1
            f.write("✓ EF 5% ≥ 5: CUMPLE\n")
        else:
            f.write(f"✗ EF 5% ≥ 5: NO CUMPLE ({metricas['ef_5']:.2f})\n")
        
        if metricas["ef_10"] >= 3:
            criterios_cumplidos += 1
            f.write("✓ EF 10% ≥ 3: CUMPLE\n")
        else:
            f.write(f"✗ EF 10% ≥ 3: NO CUMPLE ({metricas['ef_10']:.2f})\n")
        
        if metricas["bedroc_20"] >= 0.30:
            criterios_cumplidos += 1
            f.write("✓ BEDROC ≥ 0.30: CUMPLE\n")
        else:
            f.write(f"✗ BEDROC ≥ 0.30: NO CUMPLE ({metricas['bedroc_20']:.3f})\n")
        
        f.write(f"\nCriterios cumplidos: {criterios_cumplidos}/{total_criterios}\n\n")
        
        if criterios_cumplidos >= 3:
            f.write("CONCLUSIÓN: Protocolo VALIDADO para cribado virtual ACE2\n")
        else:
            f.write("CONCLUSIÓN: Protocolo SUBÓPTIMO - considerar optimización\n")
    
    # También guardar métricas en JSON
    json_file = output_path / "metricas_validacion_ace2.json"
    with open(json_file, "w") as f:
        json.dump({
            "config": CONFIG,
            "metricas": metricas,
            "fecha": datetime.now().isoformat()
        }, f, indent=2)
    
    print(f"\nReporte guardado en: {reporte_file}")
    return reporte_file


def main():
    """Función principal de validación."""
    print("\n" + "="*70)
    print("VALIDACIÓN RETROSPECTIVA GNINA - ACE2")
    print("="*70)
    print(f"Inicio: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    input_dir = Path(CONFIG["input_dir"])
    output_dir = Path(CONFIG["output_dir"])
    
    # Verificar dataset
    activos_file = input_dir / "activos.smi"
    decoys_file = input_dir / "decoys.smi"
    
    if not activos_file.exists() or not decoys_file.exists():
        print(f"\nERROR: Dataset no encontrado en {input_dir}/")
        print("Ejecuta primero:")
        print("  python obtener_activos_chembl.py")
        print("  python generar_decoys_local.py")
        return
    
    # Contar compuestos
    with open(activos_file) as f:
        n_activos = sum(1 for _ in f)
    with open(decoys_file) as f:
        n_decoys = sum(1 for _ in f)
    
    print(f"\nDataset: {n_activos} activos + {n_decoys} decoys")
    print(f"Ratio: 1:{n_decoys//max(n_activos,1)}")
    
    # Preparar directorio de salida
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Preparar ligandos
    print("\n" + "-"*40)
    print("PREPARANDO LIGANDOS...")
    print("-"*40)
    
    print("Preparando activos...")
    ligandos_activos = preparar_ligandos_pdbqt(str(activos_file), str(output_dir))
    
    print("Preparando decoys...")
    ligandos_decoys = preparar_ligandos_pdbqt(str(decoys_file), str(output_dir))
    
    print(f"\nLigandos preparados: {len(ligandos_activos)} activos, {len(ligandos_decoys)} decoys")
    
    # Crear script para ejecución en Windows con Docker
    script_windows = output_dir / "ejecutar_docking.bat"
    crear_script_windows(ligandos_activos, ligandos_decoys, str(script_windows))
    
    print(f"\nScript de ejecución creado: {script_windows}")
    print("\n" + "="*70)
    print("INSTRUCCIONES PARA EJECUCIÓN")
    print("="*70)
    print(f"""
1. Copiar carpeta {output_dir} a {CONFIG['trabajo_windows']}

2. Asegurar que el receptor está preparado:
   - Archivo: {CONFIG['receptor_file']}
   - Colocar en la misma carpeta

3. Ejecutar desde Windows (PowerShell/CMD):
   cd {CONFIG['trabajo_windows']}
   .\\ejecutar_docking.bat

4. Al completar, ejecutar cálculo de métricas:
   python calcular_metricas.py

Tiempo estimado: {(len(ligandos_activos) + len(ligandos_decoys)) * 30 // 60} minutos
(~30 seg por ligando con 2x RTX3060)
""")


def crear_script_windows(activos: List[str], decoys: List[str], output_file: str):
    """Crea script batch para Windows."""
    
    with open(output_file, "w") as f:
        f.write("@echo off\n")
        f.write("REM Script de validacion Gnina - ACE2\n")
        f.write(f"REM Generado: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        
        f.write("setlocal enabledelayedexpansion\n\n")
        
        f.write("set RECEPTOR=9FMM_receptor.pdbqt\n")
        f.write(f"set CENTER_X={CONFIG['center_x']}\n")
        f.write(f"set CENTER_Y={CONFIG['center_y']}\n")
        f.write(f"set CENTER_Z={CONFIG['center_z']}\n")
        f.write(f"set SIZE={CONFIG['size_x']}\n")
        f.write(f"set EXHAUSTIVENESS={CONFIG['exhaustiveness']}\n\n")
        
        f.write("mkdir salidas 2>nul\n")
        f.write("mkdir resultados 2>nul\n\n")
        
        f.write("echo Iniciando docking de validacion ACE2...\n")
        f.write("echo Total ligandos: %d\n\n" % (len(activos) + len(decoys)))
        
        # Activos
        f.write("echo === ACTIVOS ===\n")
        for i, lig in enumerate(activos):
            nombre = Path(lig).stem
            f.write(f"echo [{i+1}/{len(activos)}] {nombre}\n")
            f.write(f"docker run --rm --gpus all -v \"%cd%\":/data gnina/gnina gnina ")
            f.write(f"--receptor /data/%RECEPTOR% ")
            f.write(f"--ligand /data/ligandos_pdbqt/{nombre}.pdbqt ")
            f.write(f"--out /data/salidas/{nombre}_out.pdbqt ")
            f.write(f"--center_x %CENTER_X% --center_y %CENTER_Y% --center_z %CENTER_Z% ")
            f.write(f"--size_x %SIZE% --size_y %SIZE% --size_z %SIZE% ")
            f.write(f"--exhaustiveness %EXHAUSTIVENESS% ")
            f.write(f"--cnn_scoring {CONFIG['cnn_scoring']} ")
            f.write(f"--num_modes {CONFIG['num_modes']} ")
            f.write(f">> resultados/{nombre}_log.txt 2>&1\n")
        
        # Decoys
        f.write("\necho === DECOYS ===\n")
        for i, lig in enumerate(decoys):
            nombre = Path(lig).stem
            if (i + 1) % 100 == 0:
                f.write(f"echo [{i+1}/{len(decoys)}] procesando decoys...\n")
            f.write(f"docker run --rm --gpus all -v \"%cd%\":/data gnina/gnina gnina ")
            f.write(f"--receptor /data/%RECEPTOR% ")
            f.write(f"--ligand /data/ligandos_pdbqt/{nombre}.pdbqt ")
            f.write(f"--out /data/salidas/{nombre}_out.pdbqt ")
            f.write(f"--center_x %CENTER_X% --center_y %CENTER_Y% --center_z %CENTER_Z% ")
            f.write(f"--size_x %SIZE% --size_y %SIZE% --size_z %SIZE% ")
            f.write(f"--exhaustiveness %EXHAUSTIVENESS% ")
            f.write(f"--cnn_scoring {CONFIG['cnn_scoring']} ")
            f.write(f"--num_modes {CONFIG['num_modes']} ")
            f.write(f">> resultados/{nombre}_log.txt 2>&1\n")
        
        f.write("\necho.\n")
        f.write("echo Docking completado!\n")
        f.write("echo Ejecutar: python calcular_metricas.py\n")
        f.write("pause\n")


if __name__ == "__main__":
    main()
