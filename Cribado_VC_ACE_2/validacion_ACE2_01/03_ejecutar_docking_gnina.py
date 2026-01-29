#!/usr/bin/env python3
"""
03_ejecutar_docking_gnina.py
============================
Ejecuta docking molecular con Gnina vía Docker para validación ACE2.

REQUISITOS:
- Docker Desktop instalado y corriendo
- Imagen gnina/gnina descargada (docker pull gnina/gnina)
- GPU NVIDIA con drivers actualizados
- Receptor preparado: 9FMM_receptor.pdbqt
- Ligandos preparados en: resultados_validacion_ace2/ligandos_pdbqt/

PARÁMETROS DE DOCKING:
- Centro: (42.0, 7.0, 23.0) Å - sitio activo ACE2
- Caja: 26×26×26 Å
- Exhaustiveness: 64 (alto para validación)
- CNN scoring: rescore (usa scoring CNN sobre poses Vina)
- Modos: 1 (solo mejor pose por ligando)

USO:
python 03_ejecutar_docking_gnina.py

El script detecta automáticamente:
- Sistema operativo (Windows/Linux)
- Disponibilidad de GPU
- Ligandos pendientes de procesar

Autor: Script generado para TFM - Validación Gnina ACE2
Fecha: 2026-01-26
"""

import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple


# ============================================================================
# CONFIGURACIÓN
# ============================================================================

CONFIG = {
    # Receptor
    "receptor": "9FMM_receptor.pdbqt",
    
    # Centro del sitio activo (coordenadas F-MLN-4760)
    "center_x": 42.0,
    "center_y": 7.0,
    "center_z": 23.0,
    
    # Tamaño de la caja de búsqueda
    "size_x": 26,
    "size_y": 26,
    "size_z": 26,
    
    # Parámetros de búsqueda
    "exhaustiveness": 64,
    "num_modes": 1,
    "cnn_scoring": "rescore",
    
    # Directorios
    "ligandos_dir": "resultados_validacion_ace2/ligandos_pdbqt",
    "salidas_dir": "resultados_validacion_ace2/salidas_docking",
    
    # Docker
    "docker_image": "gnina/gnina",
}


# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================

def verificar_docker() -> bool:
    """Verifica que Docker esté disponible."""
    try:
        result = subprocess.run(
            ["docker", "--version"], 
            capture_output=True, 
            text=True,
            timeout=10
        )
        if result.returncode == 0:
            print(f"✓ Docker: {result.stdout.strip()}")
            return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    print("✗ Docker no disponible")
    print("  Instalar Docker Desktop: https://www.docker.com/products/docker-desktop")
    return False


def verificar_imagen_gnina() -> bool:
    """Verifica que la imagen de Gnina esté descargada."""
    try:
        result = subprocess.run(
            ["docker", "images", "-q", CONFIG["docker_image"]],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.stdout.strip():
            print(f"✓ Imagen Gnina disponible")
            return True
    except subprocess.TimeoutExpired:
        pass
    
    print(f"✗ Imagen {CONFIG['docker_image']} no encontrada")
    print(f"  Ejecutar: docker pull {CONFIG['docker_image']}")
    return False


def verificar_gpu() -> Tuple[bool, str]:
    """Verifica disponibilidad de GPU NVIDIA."""
    
    # Método 1: Verificar nvidia-smi directamente en Windows
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10
        )
        if result.returncode == 0 and result.stdout.strip():
            gpus = result.stdout.strip().split('\n')
            print(f"✓ GPU detectada: {', '.join(gpus)}")
            return True, "--gpus all"
    except (FileNotFoundError, subprocess.TimeoutExpired):
        pass
    
    # Método 2: Verificar con Docker + gnina directamente
    try:
        result = subprocess.run(
            ["docker", "run", "--rm", "--gpus", "all", 
             "gnina/gnina", "nvidia-smi", "-L"],
            capture_output=True,
            text=True,
            timeout=30
        )
        if result.returncode == 0 and "GPU" in result.stdout:
            print(f"✓ GPU disponible via Docker")
            return True, "--gpus all"
    except (subprocess.TimeoutExpired, Exception):
        pass
    
    print("⚠ GPU no disponible, usando CPU (será más lento)")
    print("  Tip: Verificar que Docker Desktop tenga GPU habilitada")
    return False, ""


def obtener_ligandos_pendientes() -> List[Path]:
    """Obtiene lista de ligandos pendientes de procesar."""
    
    ligandos_dir = Path(CONFIG["ligandos_dir"])
    salidas_dir = Path(CONFIG["salidas_dir"])
    
    if not ligandos_dir.exists():
        print(f"✗ Directorio de ligandos no encontrado: {ligandos_dir}")
        return []
    
    # Crear directorio de salidas si no existe
    salidas_dir.mkdir(parents=True, exist_ok=True)
    
    # Obtener todos los ligandos
    todos_ligandos = list(ligandos_dir.glob("*.pdbqt"))
    
    # Verificar cuáles ya están procesados
    pendientes = []
    for lig in todos_ligandos:
        salida = salidas_dir / f"{lig.stem}_out.sdf"
        if not salida.exists():
            pendientes.append(lig)
    
    return pendientes


def ejecutar_docking_ligando(ligando: Path, gpu_flag: str) -> Tuple[bool, float, str]:
    """
    Ejecuta docking para un ligando individual.
    
    Returns:
        (éxito, tiempo_segundos, mensaje)
    """
    
    salida = Path(CONFIG["salidas_dir"]) / f"{ligando.stem}_out.sdf"
    log_file = Path(CONFIG["salidas_dir"]) / f"{ligando.stem}.log"
    
    # Obtener directorio actual para montar en Docker
    work_dir = Path.cwd().absolute()
    
    # Convertir paths a formato POSIX (forward slashes) para Docker
    ligando_posix = ligando.as_posix()
    salida_posix = salida.as_posix()
    receptor_posix = CONFIG['receptor'].replace('\\', '/')
    
    # Construir comando Docker
    cmd = ["docker", "run", "--rm"]
    
    # Agregar GPU si está disponible
    if gpu_flag:
        cmd.extend(gpu_flag.split())
    
    # Montar volumen
    cmd.extend(["-v", f"{work_dir}:/data"])
    
    # Imagen y comando gnina
    cmd.extend([
        CONFIG["docker_image"],
        "gnina",
        "--receptor", f"/data/{receptor_posix}",
        "--ligand", f"/data/{ligando_posix}",
        "--out", f"/data/{salida_posix}",
        "--center_x", str(CONFIG["center_x"]),
        "--center_y", str(CONFIG["center_y"]),
        "--center_z", str(CONFIG["center_z"]),
        "--size_x", str(CONFIG["size_x"]),
        "--size_y", str(CONFIG["size_y"]),
        "--size_z", str(CONFIG["size_z"]),
        "--exhaustiveness", str(CONFIG["exhaustiveness"]),
        "--num_modes", str(CONFIG["num_modes"]),
        "--cnn_scoring", CONFIG["cnn_scoring"],
    ])
    
    # Ejecutar
    inicio = time.time()
    
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 min máximo por ligando
        )
        
        tiempo = time.time() - inicio
        
        # Guardar log
        with open(log_file, "w") as f:
            f.write(f"Comando: {' '.join(cmd)}\n")
            f.write(f"Código retorno: {result.returncode}\n")
            f.write(f"Tiempo: {tiempo:.1f} s\n")
            f.write(f"\n=== STDOUT ===\n{result.stdout}\n")
            f.write(f"\n=== STDERR ===\n{result.stderr}\n")
        
        if result.returncode == 0 and salida.exists():
            return True, tiempo, "OK"
        else:
            return False, tiempo, result.stderr[:200] if result.stderr else "Error desconocido"
            
    except subprocess.TimeoutExpired:
        return False, 300, "Timeout (>5 min)"
    except Exception as e:
        return False, 0, str(e)


def extraer_scores_sdf(sdf_path: Path) -> dict:
    """Extrae scores CNN del archivo SDF de salida."""
    
    scores = {
        "CNNscore": None,
        "CNNaffinity": None,
        "minimizedAffinity": None
    }
    
    if not sdf_path.exists():
        return scores
    
    with open(sdf_path, "r") as f:
        contenido = f.read()
    
    for linea in contenido.split('\n'):
        if "CNNscore" in linea or "> <CNNscore>" in linea:
            # Siguiente línea tiene el valor
            idx = contenido.find(linea)
            siguiente = contenido[idx:].split('\n')[1].strip()
            try:
                scores["CNNscore"] = float(siguiente)
            except:
                pass
        elif "CNNaffinity" in linea or "> <CNNaffinity>" in linea:
            idx = contenido.find(linea)
            siguiente = contenido[idx:].split('\n')[1].strip()
            try:
                scores["CNNaffinity"] = float(siguiente)
            except:
                pass
        elif "minimizedAffinity" in linea or "> <minimizedAffinity>" in linea:
            idx = contenido.find(linea)
            siguiente = contenido[idx:].split('\n')[1].strip()
            try:
                scores["minimizedAffinity"] = float(siguiente)
            except:
                pass
    
    return scores


def main():
    """Flujo principal de ejecución de docking."""
    
    print("=" * 70)
    print("DOCKING MOLECULAR CON GNINA - VALIDACIÓN ACE2")
    print("=" * 70)
    print(f"Inicio: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()
    
    # Verificaciones
    print("-" * 50)
    print("VERIFICACIÓN DE REQUISITOS")
    print("-" * 50)
    
    if not verificar_docker():
        sys.exit(1)
    
    if not verificar_imagen_gnina():
        respuesta = input("¿Descargar imagen ahora? (s/n): ")
        if respuesta.lower() == 's':
            print(f"Descargando {CONFIG['docker_image']}...")
            subprocess.run(["docker", "pull", CONFIG["docker_image"]])
        else:
            sys.exit(1)
    
    gpu_disponible, gpu_flag = verificar_gpu()
    
    # Verificar receptor
    receptor_path = Path(CONFIG["receptor"])
    if not receptor_path.exists():
        print(f"✗ Receptor no encontrado: {receptor_path}")
        print("  Ejecutar primero: python 00_preparar_receptor_ace2.py")
        sys.exit(1)
    else:
        print(f"✓ Receptor: {receptor_path}")
    
    # Obtener ligandos
    ligandos = obtener_ligandos_pendientes()
    
    if not ligandos:
        print("\n✗ No hay ligandos pendientes de procesar")
        print(f"  Verificar directorio: {CONFIG['ligandos_dir']}")
        
        # Verificar si ya están todos procesados
        salidas = list(Path(CONFIG["salidas_dir"]).glob("*_out.sdf"))
        if salidas:
            print(f"  Ya procesados: {len(salidas)} ligandos")
        
        sys.exit(0)
    
    print(f"\n✓ Ligandos pendientes: {len(ligandos)}")
    
    # Estimar tiempo
    tiempo_por_ligando = 30 if gpu_disponible else 120  # segundos
    tiempo_total_est = len(ligandos) * tiempo_por_ligando / 60
    print(f"  Tiempo estimado: {tiempo_total_est:.1f} minutos")
    
    # Confirmar ejecución
    print()
    respuesta = input("¿Iniciar docking? (s/n): ")
    if respuesta.lower() != 's':
        print("Cancelado por usuario")
        sys.exit(0)
    
    # Ejecutar docking
    print("-" * 50)
    print("EJECUTANDO DOCKING")
    print("-" * 50)
    
    resultados = []
    errores = []
    tiempo_total = 0
    
    for i, ligando in enumerate(ligandos, 1):
        print(f"[{i}/{len(ligandos)}] {ligando.stem}...", end=" ", flush=True)
        
        exito, tiempo, mensaje = ejecutar_docking_ligando(ligando, gpu_flag)
        tiempo_total += tiempo
        
        if exito:
            # Extraer scores
            salida = Path(CONFIG["salidas_dir"]) / f"{ligando.stem}_out.sdf"
            scores = extraer_scores_sdf(salida)
            
            print(f"✓ {tiempo:.1f}s | CNN: {scores['CNNscore']:.3f}" if scores['CNNscore'] else f"✓ {tiempo:.1f}s")
            
            resultados.append({
                "ligando": ligando.stem,
                "tiempo": tiempo,
                "scores": scores
            })
        else:
            print(f"✗ {mensaje}")
            errores.append({
                "ligando": ligando.stem,
                "error": mensaje
            })
    
    # Resumen final
    print()
    print("=" * 70)
    print("RESUMEN DE EJECUCIÓN")
    print("=" * 70)
    print(f"Completados: {len(resultados)}/{len(ligandos)}")
    print(f"Errores: {len(errores)}")
    print(f"Tiempo total: {tiempo_total/60:.1f} minutos")
    print(f"Promedio por ligando: {tiempo_total/max(len(ligandos),1):.1f} segundos")
    
    if errores:
        print("\nLigandos con errores:")
        for e in errores:
            print(f"  - {e['ligando']}: {e['error']}")
    
    # Guardar resultados
    resumen_path = Path(CONFIG["salidas_dir"]) / "resumen_docking.txt"
    with open(resumen_path, "w") as f:
        f.write(f"RESUMEN DOCKING ACE2 - GNINA\n")
        f.write(f"{'='*50}\n")
        f.write(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Receptor: {CONFIG['receptor']}\n")
        f.write(f"Ligandos procesados: {len(resultados)}\n")
        f.write(f"Errores: {len(errores)}\n")
        f.write(f"Tiempo total: {tiempo_total/60:.1f} min\n\n")
        
        f.write(f"PARÁMETROS:\n")
        f.write(f"  Centro: ({CONFIG['center_x']}, {CONFIG['center_y']}, {CONFIG['center_z']})\n")
        f.write(f"  Caja: {CONFIG['size_x']}×{CONFIG['size_y']}×{CONFIG['size_z']} Å\n")
        f.write(f"  Exhaustiveness: {CONFIG['exhaustiveness']}\n")
        f.write(f"  CNN scoring: {CONFIG['cnn_scoring']}\n\n")
        
        f.write(f"RESULTADOS:\n")
        f.write(f"{'Ligando':<25} {'CNNscore':>10} {'CNNaffinity':>12} {'Tiempo':>8}\n")
        f.write(f"{'-'*60}\n")
        
        for r in sorted(resultados, key=lambda x: x['scores'].get('CNNscore', 0) or 0, reverse=True):
            cnn = r['scores'].get('CNNscore', '-')
            aff = r['scores'].get('CNNaffinity', '-')
            cnn_str = f"{cnn:.4f}" if isinstance(cnn, float) else cnn
            aff_str = f"{aff:.2f}" if isinstance(aff, float) else aff
            f.write(f"{r['ligando']:<25} {cnn_str:>10} {aff_str:>12} {r['tiempo']:>7.1f}s\n")
    
    print(f"\nResultados guardados en: {resumen_path}")
    
    # Instrucciones siguientes
    print()
    print("-" * 50)
    print("SIGUIENTE PASO")
    print("-" * 50)
    print("Ejecutar cálculo de métricas:")
    print("  python calcular_metricas.py")
    
    return len(resultados), len(errores)


if __name__ == "__main__":
    exitosos, fallidos = main()
    sys.exit(0 if fallidos == 0 else 1)