#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
==============================================================================
 BENCHMARK AUTODOCK VINA - DUAL SITE (Fragment 1 + Fragment 2)
==============================================================================
 Corrección retrospectiva del benchmark original de AutoDock Vina.
 
 PROBLEMA ORIGINAL:
   El benchmark inicial usó como centro de caja las coordenadas
   (12.47, 15.46, 49.91), correspondientes al punto MEDIO (MID) entre
   los dos fragmentos idénticos del ligando A1IDX en 9FMM.cif.
   Dicho punto MID cae en espacio vacío (~24.7 Å de ambos fragmentos),
   produciendo RMSD artificialmente alto (~4.46 Å).
 
 CORRECCIÓN:
   Se implementa el enfoque dual-site (análogo al benchmark_gnina_dual_site.py),
   definiendo cajas de docking independientes centradas en cada fragmento:
     - Fragment 1 (principal): centro (16.48, 15.24, 25.54)
     - Fragment 2 (secundario): centro (8.44, 15.48, 74.27)
   Ambas cajas con dimensiones 20×20×20 Å.
   
   Además, el cálculo de RMSD emplea emparejamiento por nombre de átomo
   (name-matched Kabsch) en lugar de orden secuencial, corrigiendo el
   artefacto de ~5 Å causado por el reordenamiento atómico que introduce
   obabel al convertir PDB→PDBQT.
 
 USO:
   python benchmark_vina_dual_site.py [N_RUNS]
   python benchmark_vina_dual_site.py 5        # 5 runs por sitio
   python benchmark_vina_dual_site.py 17       # 17 runs por sitio (como Gnina)
 
 REQUISITOS:
   - Docker Desktop con imagen de Vina
   - Open Babel (obabel) en PATH
   - Python 3.8+ con numpy
   - Archivos en carpeta de trabajo:
       * 9FMM_protein.pdb (o receptor.pdbqt ya preparado)
       * ligand_exp_Fragment_1.pdb (referencia cristalográfica Frag1)
       * ligand_exp_Fragment_2.pdb (referencia cristalográfica Frag2)
       * 9FMM.cif (opcional, para verificación)
 
 HARDWARE: AMD Threadripper 3690x (48T), 128GB RAM, 2×RTX3060
 SO: Windows 10 x64
==============================================================================
"""

import subprocess
import sys
import os
import time
import json
import shutil
import re
import numpy as np
from pathlib import Path
from datetime import datetime

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

# Carpeta de trabajo
CARPETA_TRABAJO = Path(r"G:\Autodock_vina_benchmark")

# Imagen Docker de Vina (ajustar si es diferente)
DOCKER_IMAGE = "ghcr.io/metaphorme/vina:v1.2.5"

# Archivos de entrada
ARCHIVO_RECEPTOR_PDB = "9FMM_protein.pdb"
ARCHIVO_RECEPTOR_PDBQT = "receptor.pdbqt"
ARCHIVO_CIF = "9FMM.cif"

# Parámetros de docking comunes
EXHAUSTIVENESS = 64
NUM_MODES = 9
CPU_THREADS = 48
BOX_SIZE = 20  # Å (cúbico, igual para ambos sitios)

# Definición de sitios de docking
SITIOS = {
    "Fragment_1": {
        "descripcion": "Sitio principal - cadena C, A1IDX 801",
        "center_x": 16.48,
        "center_y": 15.24,
        "center_z": 25.54,
        "referencia_pdb": "ligand_exp_Fragment_1.pdb",
        "ligando_pdbqt": "ligand_Fragment_1.pdbqt",
    },
    "Fragment_2": {
        "descripcion": "Sitio secundario - cadena P, A1IDX 802",
        "center_x": 8.44,
        "center_y": 15.48,
        "center_z": 74.27,
        "referencia_pdb": "ligand_exp_Fragment_2.pdb",
        "ligando_pdbqt": "ligand_Fragment_2.pdbqt",
    },
}

# Coordenadas INCORRECTAS del benchmark original (para documentación)
CENTRO_MID_INCORRECTO = {
    "center_x": 12.470,
    "center_y": 15.461,
    "center_z": 49.907,
    "nota": "Punto MID entre fragmentos - NO corresponde a sitio de unión real"
}

# Número de runs por defecto
RUNS_POR_DEFECTO = 5


# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================

def log(mensaje, archivo_log=None):
    """Imprime y opcionalmente escribe al log."""
    marca = datetime.now().strftime("%H:%M:%S")
    linea = f"[{marca}] {mensaje}"
    print(linea)
    if archivo_log:
        with open(archivo_log, "a", encoding="utf-8") as f:
            f.write(linea + "\n")


def verificar_requisitos(archivo_log):
    """Verifica que todos los archivos y herramientas necesarios estén disponibles."""
    log("=" * 70, archivo_log)
    log("[VERIFICACIÓN] Comprobando requisitos...", archivo_log)
    log("=" * 70, archivo_log)
    
    errores = []
    
    # Verificar Open Babel
    try:
        resultado = subprocess.run(
            ["obabel", "-V"], capture_output=True, text=True, timeout=10
        )
        version_ob = resultado.stdout.strip().split("\n")[0] if resultado.stdout else "desconocida"
        log(f"  [OK] Open Babel: {version_ob}", archivo_log)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        errores.append("Open Babel (obabel) no encontrado en PATH")
        log("  [ERROR] Open Babel no encontrado", archivo_log)
    
    # Verificar Docker
    try:
        resultado = subprocess.run(
            ["docker", "info"], capture_output=True, text=True, timeout=15
        )
        if resultado.returncode == 0:
            log("  [OK] Docker disponible", archivo_log)
        else:
            errores.append("Docker no está ejecutándose")
            log("  [ERROR] Docker no responde", archivo_log)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        errores.append("Docker no encontrado")
        log("  [ERROR] Docker no encontrado", archivo_log)
    
    # Verificar receptor
    ruta_receptor_pdbqt = CARPETA_TRABAJO / ARCHIVO_RECEPTOR_PDBQT
    ruta_receptor_pdb = CARPETA_TRABAJO / ARCHIVO_RECEPTOR_PDB
    if ruta_receptor_pdbqt.exists():
        log(f"  [OK] Receptor PDBQT: {ARCHIVO_RECEPTOR_PDBQT}", archivo_log)
    elif ruta_receptor_pdb.exists():
        log(f"  [INFO] Receptor PDB encontrado, se convertirá a PDBQT", archivo_log)
    else:
        errores.append(f"Receptor no encontrado ({ARCHIVO_RECEPTOR_PDB} ni {ARCHIVO_RECEPTOR_PDBQT})")
        log(f"  [ERROR] Receptor no encontrado", archivo_log)
    
    # Verificar referencias de fragmentos
    for nombre_sitio, config in SITIOS.items():
        ruta_ref = CARPETA_TRABAJO / config["referencia_pdb"]
        if ruta_ref.exists():
            n_atomos = contar_atomos_pdb(ruta_ref)
            log(f"  [OK] Referencia {nombre_sitio}: {config['referencia_pdb']} ({n_atomos} átomos)", archivo_log)
        else:
            errores.append(f"Referencia {nombre_sitio} no encontrada: {config['referencia_pdb']}")
            log(f"  [ERROR] Referencia {nombre_sitio}: {config['referencia_pdb']} NO ENCONTRADA", archivo_log)
    
    if errores:
        log("", archivo_log)
        log("[VERIFICACIÓN] ERRORES ENCONTRADOS:", archivo_log)
        for e in errores:
            log(f"  ✗ {e}", archivo_log)
        return False
    
    log("", archivo_log)
    log("[VERIFICACIÓN] Todas las verificaciones pasaron.", archivo_log)
    return True


def contar_atomos_pdb(ruta_pdb):
    """Cuenta átomos HETATM/ATOM en un archivo PDB."""
    cuenta = 0
    with open(ruta_pdb, "r") as f:
        for linea in f:
            if linea.startswith("HETATM") or linea.startswith("ATOM"):
                cuenta += 1
    return cuenta


def extraer_coordenadas_pdb(ruta_pdb):
    """Extrae coordenadas XYZ de átomos pesados (no H) de un archivo PDB."""
    coordenadas = []
    nombres_atomos = []
    with open(ruta_pdb, "r") as f:
        for linea in f:
            if linea.startswith("HETATM") or linea.startswith("ATOM"):
                elemento = linea[76:78].strip() if len(linea) > 76 else ""
                nombre_atomo = linea[12:16].strip()
                # Excluir hidrógenos
                if elemento == "H" or (not elemento and nombre_atomo.startswith("H")):
                    continue
                try:
                    x = float(linea[30:38])
                    y = float(linea[38:46])
                    z = float(linea[46:54])
                    coordenadas.append([x, y, z])
                    nombres_atomos.append(nombre_atomo)
                except ValueError:
                    continue
    return np.array(coordenadas), nombres_atomos


def extraer_coordenadas_pdbqt(ruta_pdbqt):
    """Extrae coordenadas XYZ de átomos pesados de un archivo PDBQT (solo primer modelo)."""
    coordenadas = []
    nombres_atomos = []
    with open(ruta_pdbqt, "r") as f:
        for linea in f:
            if linea.startswith("ENDMDL"):
                break  # Solo la mejor pose (modelo 1)
            if linea.startswith("HETATM") or linea.startswith("ATOM"):
                tipo_ad = linea[77:79].strip() if len(linea) > 77 else ""
                nombre_atomo = linea[12:16].strip()
                # Excluir hidrógenos (tipo AD: HD, H)
                if tipo_ad in ("HD", "H") or (not tipo_ad and nombre_atomo.startswith("H")):
                    continue
                try:
                    x = float(linea[30:38])
                    y = float(linea[38:46])
                    z = float(linea[46:54])
                    coordenadas.append([x, y, z])
                    nombres_atomos.append(nombre_atomo)
                except ValueError:
                    continue
    return np.array(coordenadas), nombres_atomos


def calcular_rmsd_kabsch(P, Q):
    """
    Calcula RMSD usando alineamiento Kabsch (superposición óptima).
    P, Q: matrices Nx3 de coordenadas.
    Retorna RMSD en Ångströms o -1.0 si falla.
    """
    if len(P) == 0 or len(Q) == 0:
        return -1.0
    
    n = min(len(P), len(Q))
    P = P[:n]
    Q = Q[:n]
    
    # Centrar
    P_c = P - P.mean(axis=0)
    Q_c = Q - Q.mean(axis=0)
    
    # Matriz de covarianza
    H = P_c.T @ Q_c
    
    # SVD
    try:
        U, S, Vt = np.linalg.svd(H)
    except np.linalg.LinAlgError:
        return -1.0
    
    # Corrección de reflexión
    d = np.linalg.det(Vt.T @ U.T)
    signo = np.array([[1, 0, 0], [0, 1, 0], [0, 0, np.sign(d)]])
    
    # Rotación óptima
    R = Vt.T @ signo @ U.T
    
    # Aplicar rotación
    P_rot = (R @ P_c.T).T
    
    # RMSD
    diff = P_rot - Q_c
    rmsd = np.sqrt((diff ** 2).sum() / n)
    
    return rmsd


def preparar_receptor_pdbqt(archivo_log):
    """Convierte receptor PDB a PDBQT si no existe."""
    ruta_pdbqt = CARPETA_TRABAJO / ARCHIVO_RECEPTOR_PDBQT
    if ruta_pdbqt.exists():
        log(f"  Receptor PDBQT ya existe: {ARCHIVO_RECEPTOR_PDBQT}", archivo_log)
        return True
    
    ruta_pdb = CARPETA_TRABAJO / ARCHIVO_RECEPTOR_PDB
    if not ruta_pdb.exists():
        log(f"  [ERROR] Receptor PDB no encontrado", archivo_log)
        return False
    
    log(f"  Convirtiendo {ARCHIVO_RECEPTOR_PDB} → {ARCHIVO_RECEPTOR_PDBQT}...", archivo_log)
    cmd = [
        "obabel", str(ruta_pdb),
        "-O", str(ruta_pdbqt),
        "-xr",  # modo receptor
        "--partialcharge", "gasteiger"
    ]
    resultado = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    
    if ruta_pdbqt.exists() and ruta_pdbqt.stat().st_size > 0:
        log(f"  [OK] Receptor convertido", archivo_log)
        return True
    else:
        log(f"  [ERROR] Fallo en conversión: {resultado.stderr}", archivo_log)
        return False


def preparar_ligando_pdbqt(nombre_sitio, config, archivo_log):
    """Convierte ligando PDB de referencia a PDBQT para docking."""
    ruta_pdb = CARPETA_TRABAJO / config["referencia_pdb"]
    ruta_pdbqt = CARPETA_TRABAJO / config["ligando_pdbqt"]
    
    if ruta_pdbqt.exists():
        log(f"  Ligando PDBQT ya existe: {config['ligando_pdbqt']}", archivo_log)
        return True
    
    if not ruta_pdb.exists():
        log(f"  [ERROR] PDB de referencia no encontrado: {config['referencia_pdb']}", archivo_log)
        return False
    
    log(f"  Convirtiendo {config['referencia_pdb']} → {config['ligando_pdbqt']}...", archivo_log)
    cmd = [
        "obabel", str(ruta_pdb),
        "-O", str(ruta_pdbqt),
        "--partialcharge", "gasteiger",
        "-h"  # agregar hidrógenos
    ]
    resultado = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    
    if ruta_pdbqt.exists() and ruta_pdbqt.stat().st_size > 0:
        log(f"  [OK] Ligando {nombre_sitio} convertido", archivo_log)
        return True
    else:
        log(f"  [ERROR] Fallo en conversión: {resultado.stderr}", archivo_log)
        return False


def generar_config_vina(nombre_sitio, config, ruta_config):
    """Genera archivo config.txt para Vina con las coordenadas del sitio."""
    contenido = f"""# AutoDock Vina - Configuración Dual-Site
# Sitio: {nombre_sitio} ({config['descripcion']})
# Generado automáticamente por benchmark_vina_dual_site.py

receptor = /data/{ARCHIVO_RECEPTOR_PDBQT}
ligand = /data/{config['ligando_pdbqt']}

center_x = {config['center_x']:.3f}
center_y = {config['center_y']:.3f}
center_z = {config['center_z']:.3f}

size_x = {BOX_SIZE}
size_y = {BOX_SIZE}
size_z = {BOX_SIZE}

exhaustiveness = {EXHAUSTIVENESS}
num_modes = {NUM_MODES}
cpu = {CPU_THREADS}
"""
    with open(ruta_config, "w") as f:
        f.write(contenido)


def ejecutar_vina_docker(nombre_sitio, config, run_id, archivo_log):
    """Ejecuta un run de Vina en Docker y retorna ruta del output."""
    # Generar config específico para este sitio
    nombre_config = f"config_{nombre_sitio}.txt"
    ruta_config = CARPETA_TRABAJO / nombre_config
    generar_config_vina(nombre_sitio, config, ruta_config)
    
    # Nombre del archivo de salida
    nombre_salida = f"vina_out_{nombre_sitio}_{run_id:03d}.pdbqt"
    nombre_log_vina = f"vina_log_{nombre_sitio}_{run_id:03d}.txt"
    ruta_salida = CARPETA_TRABAJO / nombre_salida
    ruta_log_vina = CARPETA_TRABAJO / nombre_log_vina
    
    # Eliminar output anterior si existe
    if ruta_salida.exists():
        ruta_salida.unlink()
    
    # Construir comando Docker (sin --log, capturamos stdout)
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{CARPETA_TRABAJO}:/data",
        DOCKER_IMAGE,
        "vina",
        "--config", f"/data/{nombre_config}",
        "--out", f"/data/{nombre_salida}",
    ]
    
    log(f"  [1/3] Ejecutando Vina docking...", archivo_log)
    
    inicio = time.time()
    resultado = subprocess.run(
        cmd, capture_output=True, text=True, timeout=600
    )
    duracion = time.time() - inicio
    
    # Guardar stdout+stderr como log de Vina
    with open(ruta_log_vina, "w", encoding="utf-8") as f:
        f.write("=== STDOUT ===\n")
        f.write(resultado.stdout or "")
        f.write("\n=== STDERR ===\n")
        f.write(resultado.stderr or "")
    
    if resultado.returncode != 0:
        log(f"  [ERROR] Vina falló (código {resultado.returncode})", archivo_log)
        log(f"  stderr: {resultado.stderr[:200]}", archivo_log)
        return None, duracion, None
    
    if not ruta_salida.exists():
        log(f"  [ERROR] Archivo de salida no generado", archivo_log)
        return None, duracion, None
    
    log(f"  ✅ Docking completado en {duracion:.1f}s", archivo_log)
    
    # Extraer afinidad del output PDBQT y del log capturado
    afinidad = extraer_afinidad(ruta_salida, ruta_log_vina)
    
    return ruta_salida, duracion, afinidad


def extraer_afinidad(ruta_pdbqt, ruta_log=None):
    """Extrae la mejor afinidad (kcal/mol) del output de Vina."""
    # Intentar primero desde el output PDBQT
    try:
        with open(ruta_pdbqt, "r") as f:
            for linea in f:
                if "REMARK VINA RESULT:" in linea:
                    partes = linea.split()
                    return float(partes[3])
    except (ValueError, IndexError):
        pass
    
    # Intentar desde el log
    if ruta_log and ruta_log.exists():
        try:
            with open(ruta_log, "r") as f:
                for linea in f:
                    linea = linea.strip()
                    # Buscar línea tipo: "   1        -5.93          0          0"
                    partes = linea.split()
                    if len(partes) >= 4 and partes[0] == "1":
                        try:
                            return float(partes[1])
                        except ValueError:
                            continue
        except Exception:
            pass
    
    return None


def extraer_mejor_pose_pdbqt(ruta_pdbqt_multi):
    """Extrae la primera (mejor) pose de un PDBQT multi-modelo."""
    lineas_pose = []
    en_modelo = False
    
    with open(ruta_pdbqt_multi, "r") as f:
        for linea in f:
            if linea.startswith("MODEL"):
                en_modelo = True
                continue
            elif linea.startswith("ENDMDL"):
                break  # Solo la primera pose
            elif en_modelo or not linea.startswith("MODEL"):
                if linea.startswith("HETATM") or linea.startswith("ATOM"):
                    lineas_pose.append(linea)
    
    return lineas_pose


def emparejar_por_nombre(nombres_ref, coords_ref, nombres_dock, coords_dock):
    """
    Empareja átomos entre referencia y docked por nombre de átomo.
    Resuelve el reordenamiento que introduce obabel al convertir PDB→PDBQT.
    Retorna coordenadas emparejadas (ref, dock) y número de pares.
    """
    # Construir diccionario de referencia: nombre → coordenadas
    ref_dict = {}
    for i, nombre in enumerate(nombres_ref):
        if nombre not in ref_dict:
            ref_dict[nombre] = coords_ref[i]
    
    pareados_ref = []
    pareados_dock = []
    
    for i, nombre in enumerate(nombres_dock):
        if nombre in ref_dict:
            pareados_ref.append(ref_dict[nombre])
            pareados_dock.append(coords_dock[i])
    
    return np.array(pareados_ref), np.array(pareados_dock), len(pareados_ref)


def calcular_rmsd_run(ruta_salida_pdbqt, ruta_referencia_pdb, archivo_log):
    """
    Calcula RMSD entre la mejor pose de Vina y la referencia cristalográfica.
    Usa emparejamiento por nombre de átomo + Kabsch para evitar artefactos
    por el reordenamiento atómico que introduce obabel (PDB→PDBQT).
    """
    log(f"  [2/3] Calculando RMSD (name-matched Kabsch)...", archivo_log)
    
    # Extraer coordenadas de la referencia (PDB)
    coords_ref, nombres_ref = extraer_coordenadas_pdb(ruta_referencia_pdb)
    
    # Extraer coordenadas de la mejor pose (PDBQT)
    coords_dock, nombres_dock = extraer_coordenadas_pdbqt(ruta_salida_pdbqt)
    
    if len(coords_ref) == 0:
        log(f"  [ERROR] Sin coordenadas en referencia", archivo_log)
        return -1.0
    
    if len(coords_dock) == 0:
        log(f"  [ERROR] Sin coordenadas en pose docked", archivo_log)
        return -1.0
    
    # Emparejar átomos por nombre (corrige reordenamiento de obabel)
    ref_matched, dock_matched, n_pares = emparejar_por_nombre(
        nombres_ref, coords_ref, nombres_dock, coords_dock
    )
    
    if n_pares == 0:
        log(f"  [ERROR] No se pudieron emparejar átomos por nombre", archivo_log)
        return -1.0
    
    n_ref = len(nombres_ref)
    n_dock = len(nombres_dock)
    if n_pares < n_ref:
        log(f"  [AVISO] Emparejados {n_pares}/{n_ref} átomos pesados", archivo_log)
    
    rmsd = calcular_rmsd_kabsch(dock_matched, ref_matched)
    
    if rmsd >= 0:
        log(f"  📊 RMSD: {rmsd:.3f} Å ({n_pares} pares name-matched)", archivo_log)
    else:
        log(f"  [ERROR] Fallo en cálculo RMSD", archivo_log)
    
    return rmsd


def ejecutar_benchmark_sitio(nombre_sitio, config, n_runs, archivo_log):
    """Ejecuta el benchmark completo para un sitio."""
    log("", archivo_log)
    log("=" * 70, archivo_log)
    log(f"🎯 BENCHMARKING SITIO: {nombre_sitio}", archivo_log)
    log(f"   {config['descripcion']}", archivo_log)
    log(f"   Centro: ({config['center_x']}, {config['center_y']}, {config['center_z']})", archivo_log)
    log(f"   Caja: {BOX_SIZE}×{BOX_SIZE}×{BOX_SIZE} Å", archivo_log)
    log("=" * 70, archivo_log)
    
    # Preparar ligando PDBQT para este sitio
    if not preparar_ligando_pdbqt(nombre_sitio, config, archivo_log):
        log(f"  [ERROR FATAL] No se pudo preparar ligando para {nombre_sitio}", archivo_log)
        return None
    
    ruta_referencia = CARPETA_TRABAJO / config["referencia_pdb"]
    resultados = []
    
    for i in range(n_runs):
        log("", archivo_log)
        log(f"[{nombre_sitio}] Run {i + 1}/{n_runs}", archivo_log)
        log("-" * 50, archivo_log)
        
        hora_inicio = datetime.now().strftime("%H:%M:%S")
        log(f"  Hora inicio: {hora_inicio}", archivo_log)
        
        # Ejecutar Vina
        ruta_salida, duracion, afinidad = ejecutar_vina_docker(
            nombre_sitio, config, i + 1, archivo_log
        )
        
        hora_fin = datetime.now().strftime("%H:%M:%S")
        minutos = int(duracion // 60)
        segundos = int(duracion % 60)
        log(f"  Hora fin: {hora_fin}", archivo_log)
        log(f"  Tiempo run: {minutos}m {segundos}s ({duracion:.1f}s)", archivo_log)
        
        if ruta_salida is None:
            log(f"  [ERROR] Run fallido", archivo_log)
            resultados.append({
                "run": i + 1,
                "rmsd": -1.0,
                "afinidad": None,
                "tiempo": duracion,
                "exito": False
            })
            continue
        
        # Calcular RMSD
        rmsd = calcular_rmsd_run(ruta_salida, ruta_referencia, archivo_log)
        
        if afinidad is not None:
            log(f"  Afinidad: {afinidad:.3f} kcal/mol", archivo_log)
        
        resultados.append({
            "run": i + 1,
            "rmsd": rmsd,
            "afinidad": afinidad,
            "tiempo": duracion,
            "exito": rmsd > 0
        })
    
    return resultados


def calcular_estadisticas(resultados, nombre_sitio, archivo_log):
    """Calcula y reporta estadísticas para un sitio."""
    exitosos = [r for r in resultados if r["exito"]]
    n_exitosos = len(exitosos)
    n_total = len(resultados)
    
    log("", archivo_log)
    log(f"📊 ESTADÍSTICAS {nombre_sitio} (n={n_exitosos}/{n_total}):", archivo_log)
    log("-" * 50, archivo_log)
    
    if n_exitosos == 0:
        log(f"  ⚠️ Sin runs exitosos para {nombre_sitio}", archivo_log)
        return None
    
    rmsds = [r["rmsd"] for r in exitosos]
    tiempos = [r["tiempo"] for r in exitosos]
    afinidades = [r["afinidad"] for r in exitosos if r["afinidad"] is not None]
    
    stats = {
        "nombre": nombre_sitio,
        "n_total": n_total,
        "n_exitosos": n_exitosos,
        "rmsd_promedio": np.mean(rmsds),
        "rmsd_std": np.std(rmsds),
        "rmsd_min": np.min(rmsds),
        "rmsd_max": np.max(rmsds),
        "rmsd_mediana": np.median(rmsds),
        "success_2A": sum(1 for r in rmsds if r <= 2.0),
        "success_3A": sum(1 for r in rmsds if r <= 3.0),
        "tiempo_promedio": np.mean(tiempos),
        "tiempo_total": sum(tiempos),
        "afinidad_promedio": np.mean(afinidades) if afinidades else None,
        "afinidad_std": np.std(afinidades) if afinidades else None,
    }
    
    log(f"  RMSD promedio:    {stats['rmsd_promedio']:.3f} ± {stats['rmsd_std']:.3f} Å", archivo_log)
    log(f"  RMSD rango:       {stats['rmsd_min']:.3f} - {stats['rmsd_max']:.3f} Å", archivo_log)
    log(f"  RMSD mediana:     {stats['rmsd_mediana']:.3f} Å", archivo_log)
    log(f"  Success ≤2.0 Å:   {stats['success_2A']}/{n_exitosos} ({100*stats['success_2A']/n_exitosos:.1f}%)", archivo_log)
    log(f"  Success ≤3.0 Å:   {stats['success_3A']}/{n_exitosos} ({100*stats['success_3A']/n_exitosos:.1f}%)", archivo_log)
    
    if stats["afinidad_promedio"] is not None:
        log(f"  Afinidad prom.:   {stats['afinidad_promedio']:.3f} ± {stats['afinidad_std']:.3f} kcal/mol", archivo_log)
    
    minutos_total = int(stats["tiempo_total"] // 60)
    segundos_total = int(stats["tiempo_total"] % 60)
    log(f"  Tiempo promedio:  {stats['tiempo_promedio']:.1f} s/run", archivo_log)
    log(f"  Tiempo total:     {minutos_total}m {segundos_total}s ({stats['tiempo_total']:.1f}s)", archivo_log)
    
    return stats


def tabla_resultados_por_run(resultados, nombre_sitio, archivo_log):
    """Imprime tabla detallada de resultados por run."""
    log("", archivo_log)
    log(f"RESULTADOS POR RUN - {nombre_sitio}:", archivo_log)
    log("-" * 70, archivo_log)
    log(f"{'Run':<6}{'Tiempo(s)':<12}{'RMSD(Å)':<12}{'Afinidad(kcal/mol)':<22}{'Estado'}", archivo_log)
    log("-" * 70, archivo_log)
    
    for r in resultados:
        rmsd_str = f"{r['rmsd']:.3f}" if r["rmsd"] > 0 else "FALLO"
        afin_str = f"{r['afinidad']:.3f}" if r["afinidad"] is not None else "N/A"
        estado = "✅" if r["exito"] else "❌"
        log(f"{r['run']:<6}{r['tiempo']:<12.1f}{rmsd_str:<12}{afin_str:<22}{estado}", archivo_log)
    
    log("-" * 70, archivo_log)


def comparacion_dual_site(stats_sitios, archivo_log):
    """Compara resultados entre ambos sitios."""
    log("", archivo_log)
    log("=" * 70, archivo_log)
    log("📊 COMPARACIÓN DUAL-SITE", archivo_log)
    log("=" * 70, archivo_log)
    
    sitios_validos = {k: v for k, v in stats_sitios.items() if v is not None}
    
    if len(sitios_validos) == 0:
        log("  ⚠️ Sin resultados válidos para comparar", archivo_log)
        return
    
    # Tabla comparativa
    log("", archivo_log)
    log(f"{'Métrica':<25}{'Fragment_1':<18}{'Fragment_2':<18}{'Diferencia'}", archivo_log)
    log("-" * 75, archivo_log)
    
    if len(sitios_validos) == 2:
        s1 = sitios_validos.get("Fragment_1")
        s2 = sitios_validos.get("Fragment_2")
        
        if s1 and s2:
            log(f"{'RMSD promedio':<25}{s1['rmsd_promedio']:<18.3f}{s2['rmsd_promedio']:<18.3f}{abs(s1['rmsd_promedio']-s2['rmsd_promedio']):.3f} Å", archivo_log)
            log(f"{'RMSD desv. est.':<25}{s1['rmsd_std']:<18.3f}{s2['rmsd_std']:<18.3f}{abs(s1['rmsd_std']-s2['rmsd_std']):.3f} Å", archivo_log)
            log(f"{'Success ≤2.0 Å':<25}{s1['success_2A']}/{s1['n_exitosos']:<15}{s2['success_2A']}/{s2['n_exitosos']:<15}", archivo_log)
            log(f"{'Success ≤3.0 Å':<25}{s1['success_3A']}/{s1['n_exitosos']:<15}{s2['success_3A']}/{s2['n_exitosos']:<15}", archivo_log)
            log(f"{'Tiempo prom. (s)':<25}{s1['tiempo_promedio']:<18.1f}{s2['tiempo_promedio']:<18.1f}", archivo_log)
            
            # Mejor sitio
            mejor = "Fragment_1" if s1["rmsd_promedio"] < s2["rmsd_promedio"] else "Fragment_2"
            mejor_rmsd = min(s1["rmsd_promedio"], s2["rmsd_promedio"])
            diff = abs(s1["rmsd_promedio"] - s2["rmsd_promedio"])
            
            log("", archivo_log)
            log(f"🏆 MEJOR SITIO: {mejor}", archivo_log)
            log(f"   RMSD promedio: {mejor_rmsd:.3f} Å", archivo_log)
            log(f"   Diferencia: {diff:.3f} Å mejor que el otro sitio", archivo_log)
    else:
        for nombre, stats in sitios_validos.items():
            log(f"\n  {nombre}:", archivo_log)
            log(f"    RMSD: {stats['rmsd_promedio']:.3f} ± {stats['rmsd_std']:.3f} Å", archivo_log)
            log(f"    Success ≤2.0 Å: {stats['success_2A']}/{stats['n_exitosos']}", archivo_log)
    
    # Comparación con benchmark original (punto MID + RMSD sin name-match)
    log("", archivo_log)
    log("=" * 70, archivo_log)
    log("📋 COMPARACIÓN CON BENCHMARK ORIGINAL", archivo_log)
    log("=" * 70, archivo_log)
    log(f"  Benchmark original tenía DOS errores:", archivo_log)
    log(f"    1) Centro de caja en punto MID: ({CENTRO_MID_INCORRECTO['center_x']}, "
        f"{CENTRO_MID_INCORRECTO['center_y']}, {CENTRO_MID_INCORRECTO['center_z']})", archivo_log)
    log(f"    2) RMSD por orden secuencial (sin name-match) → artefacto ~5 Å", archivo_log)
    log(f"  RMSD original reportado:  4.446 ± 0.140 Å (exh=64, n=5)", archivo_log)
    log("", archivo_log)
    
    for nombre, stats in sitios_validos.items():
        mejora = 4.446 - stats["rmsd_promedio"]
        mejora_pct = (mejora / 4.446) * 100
        log(f"  {nombre} corregido:    {stats['rmsd_promedio']:.3f} Å "
            f"(mejora: {mejora:.3f} Å / {mejora_pct:.1f}%)", archivo_log)
    
    # Calificación
    log("", archivo_log)
    for nombre, stats in sitios_validos.items():
        rmsd = stats["rmsd_promedio"]
        if rmsd <= 1.0:
            calificacion = "⭐ EXCEPCIONAL"
        elif rmsd <= 2.0:
            calificacion = "✅ EXCELENTE"
        elif rmsd <= 3.0:
            calificacion = "🟡 ACEPTABLE"
        else:
            calificacion = "❌ INSUFICIENTE"
        log(f"  {nombre}: {rmsd:.3f} Å → {calificacion}", archivo_log)


def guardar_resultados_json(todos_resultados, stats_sitios, n_runs, archivo_json):
    """Guarda todos los resultados en formato JSON para análisis posterior."""
    datos = {
        "metadata": {
            "script": "benchmark_vina_dual_site.py",
            "fecha": datetime.now().isoformat(),
            "tipo": "corrección retrospectiva dual-site",
            "programa": "AutoDock Vina",
            "exhaustiveness": EXHAUSTIVENESS,
            "num_modes": NUM_MODES,
            "cpu_threads": CPU_THREADS,
            "box_size": BOX_SIZE,
            "n_runs_por_sitio": n_runs,
            "centro_incorrecto_original": CENTRO_MID_INCORRECTO,
        },
        "sitios": {},
    }
    
    for nombre_sitio, resultados in todos_resultados.items():
        stats = stats_sitios.get(nombre_sitio)
        datos["sitios"][nombre_sitio] = {
            "config": {
                "center_x": SITIOS[nombre_sitio]["center_x"],
                "center_y": SITIOS[nombre_sitio]["center_y"],
                "center_z": SITIOS[nombre_sitio]["center_z"],
            },
            "resultados": resultados,
            "estadisticas": {
                k: float(v) if isinstance(v, (np.floating, np.integer)) else v
                for k, v in (stats or {}).items()
            } if stats else None,
        }
    
    with open(archivo_json, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=2, ensure_ascii=False, default=str)


# ============================================================================
# FUNCIÓN PRINCIPAL
# ============================================================================

def main():
    # Número de runs desde CLI
    n_runs = RUNS_POR_DEFECTO
    if len(sys.argv) > 1:
        try:
            n_runs = int(sys.argv[1])
        except ValueError:
            print(f"Uso: python {sys.argv[0]} [N_RUNS]")
            print(f"  N_RUNS: número de runs por sitio (default: {RUNS_POR_DEFECTO})")
            sys.exit(1)
    
    # Crear carpeta de resultados
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    carpeta_resultados = CARPETA_TRABAJO / f"resultados_dual_site_{timestamp}"
    carpeta_resultados.mkdir(exist_ok=True)
    
    # Archivo de log consolidado
    archivo_log = carpeta_resultados / "benchmark_dual_site.log"
    
    # Banner
    log("=" * 70, archivo_log)
    log("🧬 AUTODOCK VINA - BENCHMARK DUAL-SITE (CORRECCIÓN RETROSPECTIVA)", archivo_log)
    log("=" * 70, archivo_log)
    log(f"  Target:         9FMM (ACE2 + F-MLN-4760)", archivo_log)
    log(f"  Programa:       AutoDock Vina (Docker: {DOCKER_IMAGE})", archivo_log)
    log(f"  Exhaustiveness: {EXHAUSTIVENESS}", archivo_log)
    log(f"  CPU threads:    {CPU_THREADS}", archivo_log)
    log(f"  Box size:       {BOX_SIZE}×{BOX_SIZE}×{BOX_SIZE} Å", archivo_log)
    log(f"  Runs por sitio: {n_runs}", archivo_log)
    log(f"  Total runs:     {n_runs * len(SITIOS)}", archivo_log)
    log(f"  Fecha:          {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", archivo_log)
    log(f"  Resultados en:  {carpeta_resultados}", archivo_log)
    log("", archivo_log)
    log("  CORRECCIÓN 1: Centro de caja → de punto MID a fragmentos individuales", archivo_log)
    log(f"    MID (incorrecto):  ({CENTRO_MID_INCORRECTO['center_x']}, "
        f"{CENTRO_MID_INCORRECTO['center_y']}, {CENTRO_MID_INCORRECTO['center_z']})", archivo_log)
    log(f"    Fragment 1:        ({SITIOS['Fragment_1']['center_x']}, "
        f"{SITIOS['Fragment_1']['center_y']}, {SITIOS['Fragment_1']['center_z']})", archivo_log)
    log(f"    Fragment 2:        ({SITIOS['Fragment_2']['center_x']}, "
        f"{SITIOS['Fragment_2']['center_y']}, {SITIOS['Fragment_2']['center_z']})", archivo_log)
    log("  CORRECCIÓN 2: RMSD → name-matched Kabsch (corrige reorden obabel)", archivo_log)
    log("=" * 70, archivo_log)
    
    # Verificaciones
    if not verificar_requisitos(archivo_log):
        log("\n[FATAL] Verificaciones fallidas. Corregir antes de continuar.", archivo_log)
        sys.exit(1)
    
    # Preparar receptor
    log("", archivo_log)
    log("[PREPARACIÓN] Preparando receptor...", archivo_log)
    if not preparar_receptor_pdbqt(archivo_log):
        log("[FATAL] No se pudo preparar receptor.", archivo_log)
        sys.exit(1)
    
    # Ejecutar benchmarks para cada sitio
    inicio_global = time.time()
    todos_resultados = {}
    stats_sitios = {}
    
    for nombre_sitio, config in SITIOS.items():
        resultados = ejecutar_benchmark_sitio(nombre_sitio, config, n_runs, archivo_log)
        
        if resultados is None:
            log(f"\n[ERROR] Benchmark {nombre_sitio} falló completamente", archivo_log)
            todos_resultados[nombre_sitio] = []
            stats_sitios[nombre_sitio] = None
            continue
        
        todos_resultados[nombre_sitio] = resultados
        
        # Tabla detallada
        tabla_resultados_por_run(resultados, nombre_sitio, archivo_log)
        
        # Estadísticas
        stats = calcular_estadisticas(resultados, nombre_sitio, archivo_log)
        stats_sitios[nombre_sitio] = stats
    
    duracion_global = time.time() - inicio_global
    
    # Comparación dual-site
    comparacion_dual_site(stats_sitios, archivo_log)
    
    # Tiempo total
    min_total = int(duracion_global // 60)
    seg_total = int(duracion_global % 60)
    log("", archivo_log)
    log(f"⏱️  TIEMPO TOTAL BENCHMARK: {min_total}m {seg_total}s ({duracion_global:.1f}s)", archivo_log)
    
    # Guardar JSON
    archivo_json = carpeta_resultados / "resultados_dual_site.json"
    guardar_resultados_json(todos_resultados, stats_sitios, n_runs, archivo_json)
    log(f"\n📁 Resultados guardados en:", archivo_log)
    log(f"   Log:  {archivo_log}", archivo_log)
    log(f"   JSON: {archivo_json}", archivo_log)
    
    # Copiar este script como referencia
    try:
        shutil.copy2(__file__, carpeta_resultados / "script_usado.py")
    except Exception:
        pass
    
    log("\n✅ Benchmark dual-site completado.", archivo_log)


if __name__ == "__main__":
    main()