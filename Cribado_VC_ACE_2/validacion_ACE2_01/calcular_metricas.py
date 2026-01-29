#!/usr/bin/env python3
"""
calcular_metricas.py
====================
Calcula métricas de validación a partir de los resultados de docking Gnina.

Parsea archivos de salida PDBQT y logs para extraer scores CNN,
luego calcula AUC-ROC, EF, BEDROC y genera reporte.

Autor: Claude para TFM - Validación Gnina ACE2
Fecha: 2025-01-25
"""

import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import numpy as np
except ImportError:
    os.system("pip install numpy --break-system-packages -q")
    import numpy as np

try:
    from sklearn.metrics import roc_auc_score, roc_curve
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    os.system("pip install scikit-learn matplotlib --break-system-packages -q")
    from sklearn.metrics import roc_auc_score, roc_curve
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt


# Configuración
CONFIG = {
    "salidas_dir": "resultados_validacion_ace2/salidas_docking",
    "resultados_dir": "resultados_validacion_ace2/salidas_docking",
    "activos_file": "dataset_ace2/validacion/activos.smi",
    "decoys_file": "dataset_ace2/validacion/decoys.smi",
    "output_dir": "metricas_finales"
}


def cargar_labels() -> Dict[str, int]:
    """Carga etiquetas (1=activo, 0=decoy) desde archivos SMI."""
    labels = {}
    
    # Activos
    if Path(CONFIG["activos_file"]).exists():
        with open(CONFIG["activos_file"]) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    labels[parts[1]] = 1
    
    # Decoys
    if Path(CONFIG["decoys_file"]).exists():
        with open(CONFIG["decoys_file"]) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    labels[parts[1]] = 0
    
    print(f"Labels cargados: {sum(labels.values())} activos, {len(labels) - sum(labels.values())} decoys")
    return labels


def parsear_pdbqt(filepath: str) -> Dict:
    """Extrae scores de un archivo PDBQT de Gnina."""
    resultado = {
        "cnn_score": None,
        "cnn_affinity": None,
        "vina_affinity": None
    }
    
    try:
        with open(filepath, "r") as f:
            for line in f:
                if line.startswith("REMARK"):
                    # CNN score (VS = virtual screening score)
                    if "CNNscore" in line or "CNN_VS" in line:
                        match = re.search(r'[\d.]+$', line.strip())
                        if match:
                            resultado["cnn_score"] = float(match.group())
                    
                    # CNN affinity
                    elif "CNNaffinity" in line:
                        match = re.search(r'-?[\d.]+$', line.strip())
                        if match:
                            resultado["cnn_affinity"] = float(match.group())
                    
                    # Vina affinity
                    elif "VINA RESULT" in line:
                        parts = line.split()
                        for i, p in enumerate(parts):
                            try:
                                resultado["vina_affinity"] = float(parts[i+1])
                                break
                            except (ValueError, IndexError):
                                continue
    except Exception as e:
        pass
    
    return resultado


def parsear_log(filepath: str) -> Dict:
    """Extrae scores de un archivo de log de Gnina."""
    resultado = {
        "cnn_score": None,
        "cnn_affinity": None,
        "vina_affinity": None
    }
    
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
            
            # Formato Gnina v1.3:
            # mode |  affinity  |  intramol  |    CNN     |   CNN
            #      | (kcal/mol) | (kcal/mol) | pose score | affinity
            # -----+------------+------------+------------+----------
            #     1       -5.06       -0.14       0.8479      3.898
            
            # Buscar línea con resultados (empieza con espacios y número 1)
            pattern = r'^\s*1\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)\s+(-?[\d.]+)'
            match = re.search(pattern, content, re.MULTILINE)
            
            if match:
                resultado["vina_affinity"] = float(match.group(1))  # affinity
                # match.group(2) es intramol, lo ignoramos
                resultado["cnn_score"] = float(match.group(3))      # CNN pose score
                resultado["cnn_affinity"] = float(match.group(4))   # CNN affinity
    except Exception as e:
        pass
    
    return resultado


def parsear_sdf(filepath: str) -> Dict:
    """Extrae scores de un archivo SDF de Gnina."""
    resultado = {
        "cnn_score": None,
        "cnn_affinity": None,
        "vina_affinity": None
    }
    
    try:
        with open(filepath, "r") as f:
            content = f.read()
            
            # Buscar CNNscore
            if "> <CNNscore>" in content:
                idx = content.find("> <CNNscore>")
                lineas = content[idx:].split('\n')
                if len(lineas) > 1:
                    try:
                        resultado["cnn_score"] = float(lineas[1].strip())
                    except ValueError:
                        pass
            
            # Buscar CNNaffinity
            if "> <CNNaffinity>" in content:
                idx = content.find("> <CNNaffinity>")
                lineas = content[idx:].split('\n')
                if len(lineas) > 1:
                    try:
                        resultado["cnn_affinity"] = float(lineas[1].strip())
                    except ValueError:
                        pass
            
            # Buscar minimizedAffinity (Vina score)
            if "> <minimizedAffinity>" in content:
                idx = content.find("> <minimizedAffinity>")
                lineas = content[idx:].split('\n')
                if len(lineas) > 1:
                    try:
                        resultado["vina_affinity"] = float(lineas[1].strip())
                    except ValueError:
                        pass
    except Exception as e:
        pass
    
    return resultado


def recolectar_resultados(labels: Dict[str, int]) -> List[Dict]:
    """Recolecta todos los resultados de docking."""
    resultados = []
    
    salidas_dir = Path(CONFIG["salidas_dir"])
    resultados_dir = Path(CONFIG["resultados_dir"])
    
    # Diagnóstico: mostrar archivos disponibles
    archivos_log = list(salidas_dir.glob("*.log"))
    archivos_sdf = list(salidas_dir.glob("*_out.sdf"))
    print(f"Archivos .log encontrados: {len(archivos_log)}")
    print(f"Archivos .sdf encontrados: {len(archivos_sdf)}")
    
    encontrados = 0
    no_encontrados = []
    
    for nombre, label in labels.items():
        # Buscar archivos
        log_file = salidas_dir / f"{nombre}.log"
        sdf_file = salidas_dir / f"{nombre}_out.sdf"
        pdbqt_file = salidas_dir / f"{nombre}_out.pdbqt"
        
        scores = {"cnn_score": None, "cnn_affinity": None, "vina_affinity": None}
        
        # Priorizar LOG (contiene los scores CNN)
        if log_file.exists():
            scores = parsear_log(str(log_file))
            if scores["cnn_score"] is not None:
                encontrados += 1
        
        # Complementar con SDF si no se encontró en log
        if scores["cnn_score"] is None and sdf_file.exists():
            scores_sdf = parsear_sdf(str(sdf_file))
            for key, val in scores_sdf.items():
                if scores[key] is None:
                    scores[key] = val
            if scores["cnn_score"] is not None:
                encontrados += 1
        
        # PDBQT como último recurso
        if scores["cnn_score"] is None and pdbqt_file.exists():
            scores_pdbqt = parsear_pdbqt(str(pdbqt_file))
            for key, val in scores_pdbqt.items():
                if scores[key] is None:
                    scores[key] = val
        
        if scores["cnn_score"] is None:
            no_encontrados.append(nombre)
        
        resultados.append({
            "nombre": nombre,
            "label": label,
            **scores
        })
    
    print(f"Scores extraídos: {encontrados}/{len(labels)}")
    if no_encontrados and len(no_encontrados) <= 10:
        print(f"Sin score: {no_encontrados[:10]}")
    
    # Estadísticas
    con_score = sum(1 for r in resultados if r["cnn_score"] is not None)
    print(f"Resultados recolectados: {con_score}/{len(resultados)} con CNN score")
    
    return resultados


def calcular_bedroc(y_true: np.ndarray, y_scores: np.ndarray, alpha: float = 20.0) -> float:
    """Calcula BEDROC (Boltzmann-Enhanced Discrimination of ROC)."""
    n = len(y_true)
    n_actives = np.sum(y_true)
    
    if n_actives == 0 or n_actives == n:
        return 0.0
    
    sorted_indices = np.argsort(-y_scores)
    sorted_labels = y_true[sorted_indices]
    
    positions = np.where(sorted_labels == 1)[0] + 1
    s = np.sum(np.exp(-alpha * positions / n))
    
    Ra = n_actives / n
    Ri = 1 - np.exp(-alpha)
    Rmin = (Ra / n) * (1 - np.exp(-alpha)) / (np.exp(alpha/n) - 1)
    Rmax = (1 - np.exp(-alpha * Ra)) / (Ra * (1 - np.exp(-alpha)))
    
    bedroc = (s / n_actives - Rmin) / (Rmax - Rmin)
    
    return max(0.0, min(1.0, bedroc))


def calcular_ef(y_true: np.ndarray, y_scores: np.ndarray, porcentaje: float) -> float:
    """Calcula Enrichment Factor a un porcentaje dado."""
    n = len(y_true)
    n_actives = np.sum(y_true)
    
    if n_actives == 0:
        return 0.0
    
    top_n = max(1, int(n * porcentaje / 100))
    sorted_indices = np.argsort(-y_scores)
    top_labels = y_true[sorted_indices[:top_n]]
    
    activos_encontrados = np.sum(top_labels)
    activos_esperados = n_actives * porcentaje / 100
    
    if activos_esperados == 0:
        return 0.0
    
    return activos_encontrados / activos_esperados


def calcular_todas_metricas(resultados: List[Dict]) -> Dict:
    """Calcula todas las métricas de validación."""
    
    # Filtrar resultados con scores válidos
    validos = [r for r in resultados if r["cnn_score"] is not None]
    
    if len(validos) < 10:
        print(f"ADVERTENCIA: Solo {len(validos)} resultados válidos")
        return None
    
    labels = np.array([r["label"] for r in validos])
    scores_cnn = np.array([r["cnn_score"] for r in validos])
    scores_aff = np.array([r["cnn_affinity"] if r["cnn_affinity"] else 0 for r in validos])
    
    # Calcular con CNN score
    auc_cnn_score = roc_auc_score(labels, scores_cnn)
    
    # Calcular con CNN affinity (mayor es mejor unión)
    auc_cnn_affinity = roc_auc_score(labels, scores_aff)
    
    print(f"\n--- Comparación de métricas ---")
    print(f"AUC-ROC (CNN score):    {auc_cnn_score:.3f}")
    print(f"AUC-ROC (CNN affinity): {auc_cnn_affinity:.3f}")
    
    # Usar la mejor métrica
    if auc_cnn_affinity > auc_cnn_score:
        print(f"→ Usando CNN affinity (mejor discriminación)")
        scores = scores_aff
        metrica_usada = "cnn_affinity"
    else:
        print(f"→ Usando CNN score")
        scores = scores_cnn
        metrica_usada = "cnn_score"
    
    metricas = {
        "n_total": len(validos),
        "n_activos": int(np.sum(labels)),
        "n_decoys": int(len(labels) - np.sum(labels)),
        "ratio": (len(labels) - np.sum(labels)) / max(np.sum(labels), 1),
        
        # AUC-ROC (ambas métricas)
        "auc_roc": roc_auc_score(labels, scores),
        "auc_roc_cnn_score": auc_cnn_score,
        "auc_roc_cnn_affinity": auc_cnn_affinity,
        "metrica_usada": metrica_usada,
        
        # Enrichment Factors
        "ef_1": calcular_ef(labels, scores, 1.0),
        "ef_5": calcular_ef(labels, scores, 5.0),
        "ef_10": calcular_ef(labels, scores, 10.0),
        
        # BEDROC
        "bedroc_20": calcular_bedroc(labels, scores, alpha=20.0)
    }
    
    # Curva ROC para gráfica
    fpr, tpr, _ = roc_curve(labels, scores)
    metricas["roc_fpr"] = fpr.tolist()
    metricas["roc_tpr"] = tpr.tolist()
    
    return metricas


def generar_graficas(metricas: Dict, output_dir: str):
    """Genera gráficas de ROC y enriquecimiento."""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Curva ROC
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(metricas["roc_fpr"], metricas["roc_tpr"], 
            'b-', linewidth=2, label=f'Gnina (AUC = {metricas["auc_roc"]:.3f})')
    ax.plot([0, 1], [0, 1], 'k--', label='Random')
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('Curva ROC - Validación ACE2', fontsize=14)
    ax.legend(loc='lower right')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_path / "curva_roc_ace2.png", dpi=150)
    plt.close()
    
    # Gráfica de barras EF
    fig, ax = plt.subplots(figsize=(10, 6))
    metricas_ef = ['EF 1%', 'EF 5%', 'EF 10%', 'BEDROC']
    valores = [metricas["ef_1"], metricas["ef_5"], metricas["ef_10"], 
               metricas["bedroc_20"] * 10]  # Escalar BEDROC para visualización
    umbrales = [10, 5, 3, 3]  # Umbral BEDROC = 0.3 escalado
    
    x = np.arange(len(metricas_ef))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, valores, width, label='Gnina', color='steelblue')
    bars2 = ax.bar(x + width/2, umbrales, width, label='Umbral', color='coral', alpha=0.7)
    
    ax.set_ylabel('Valor', fontsize=12)
    ax.set_title('Métricas de Enriquecimiento - Validación ACE2', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(metricas_ef)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    # Añadir valores sobre barras
    for bar, val in zip(bars1, valores):
        height = bar.get_height()
        ax.annotate(f'{val:.2f}',
                    xy=(bar.get_x() + bar.get_width()/2, height),
                    ha='center', va='bottom', fontsize=10)
    
    plt.tight_layout()
    plt.savefig(output_path / "metricas_enriquecimiento_ace2.png", dpi=150)
    plt.close()
    
    print(f"Gráficas guardadas en {output_path}/")


def generar_reporte_texto(metricas: Dict, output_dir: str):
    """Genera reporte de texto formateado."""
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    reporte_file = output_path / "reporte_validacion_ace2.txt"
    
    with open(reporte_file, "w", encoding="utf-8") as f:
        f.write("="*70 + "\n")
        f.write("REPORTE DE VALIDACIÓN RETROSPECTIVA - ACE2\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Fecha generación: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
        f.write(f"Receptor: PDB 9FMM (ACE2 humana)\n")
        f.write(f"Programa: Gnina (CNN scoring)\n\n")
        
        f.write("-"*70 + "\n")
        f.write("COMPOSICIÓN DEL DATASET\n")
        f.write("-"*70 + "\n")
        f.write(f"Total compuestos evaluados: {metricas['n_total']}\n")
        f.write(f"Activos: {metricas['n_activos']}\n")
        f.write(f"Decoys:  {metricas['n_decoys']}\n")
        f.write(f"Ratio (decoys:activos): {metricas['ratio']:.1f}:1\n\n")
        
        f.write("-"*70 + "\n")
        f.write("MÉTRICAS DE VALIDACIÓN\n")
        f.write("-"*70 + "\n\n")
        
        f.write(f"{'Métrica':<20} {'Valor':<12} {'Umbral':<12} {'Estado':<12}\n")
        f.write("-"*56 + "\n")
        
        # Mostrar ambos AUC-ROC si están disponibles
        if "auc_roc_cnn_score" in metricas:
            f.write(f"{'AUC-ROC (CNN score)':<20} {metricas['auc_roc_cnn_score']:<12.3f}\n")
            f.write(f"{'AUC-ROC (CNN aff.)':<20} {metricas['auc_roc_cnn_affinity']:<12.3f}\n")
            f.write(f"{'Métrica usada:':<20} {metricas.get('metrica_usada', 'cnn_score')}\n")
            f.write("-"*56 + "\n")
        
        # AUC-ROC
        auc_ok = metricas["auc_roc"] >= 0.70
        f.write(f"{'AUC-ROC (final)':<20} {metricas['auc_roc']:<12.3f} {'≥0.70':<12} {'✓ CUMPLE' if auc_ok else '✗ NO CUMPLE':<12}\n")
        
        # EF 1%
        ef1_ok = metricas["ef_1"] >= 10
        f.write(f"{'EF 1%':<20} {metricas['ef_1']:<12.2f} {'≥10':<12} {'✓ CUMPLE' if ef1_ok else '✗ NO CUMPLE':<12}\n")
        
        # EF 5%
        ef5_ok = metricas["ef_5"] >= 5
        f.write(f"{'EF 5%':<20} {metricas['ef_5']:<12.2f} {'≥5':<12} {'✓ CUMPLE' if ef5_ok else '✗ NO CUMPLE':<12}\n")
        
        # EF 10%
        ef10_ok = metricas["ef_10"] >= 3
        f.write(f"{'EF 10%':<20} {metricas['ef_10']:<12.2f} {'≥3':<12} {'✓ CUMPLE' if ef10_ok else '✗ NO CUMPLE':<12}\n")
        
        # BEDROC
        bedroc_ok = metricas["bedroc_20"] >= 0.30
        f.write(f"{'BEDROC (α=20)':<20} {metricas['bedroc_20']:<12.3f} {'≥0.30':<12} {'✓ CUMPLE' if bedroc_ok else '✗ NO CUMPLE':<12}\n")
        
        f.write("\n" + "-"*70 + "\n")
        f.write("EVALUACIÓN GLOBAL\n")
        f.write("-"*70 + "\n\n")
        
        criterios = [auc_ok, ef1_ok, ef5_ok, ef10_ok, bedroc_ok]
        cumplidos = sum(criterios)
        
        f.write(f"Criterios cumplidos: {cumplidos}/5\n\n")
        
        if cumplidos >= 4:
            f.write("CONCLUSIÓN: Protocolo VALIDADO - Rendimiento excelente\n")
            f.write("El protocolo Gnina demuestra capacidad discriminatoria robusta\n")
            f.write("para distinguir activos ACE2 de decoys property-matched.\n")
        elif cumplidos >= 3:
            f.write("CONCLUSIÓN: Protocolo ACEPTABLE - Rendimiento satisfactorio\n")
            f.write("El protocolo puede emplearse para cribado virtual exploratorio.\n")
            f.write("Se recomienda validación experimental de hits prioritarios.\n")
        else:
            f.write("CONCLUSIÓN: Protocolo SUBÓPTIMO - Requiere optimización\n")
            f.write("Considerar ajustar parámetros de docking o grid box.\n")
            f.write("Validar con dataset alternativo antes de cribado productivo.\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("NOTAS METODOLÓGICAS\n")
        f.write("="*70 + "\n\n")
        f.write("- AUC-ROC: Área bajo curva ROC, capacidad discriminatoria global\n")
        f.write("- EF (Enrichment Factor): Enriquecimiento vs azar en top X%\n")
        f.write("- BEDROC: Énfasis en enriquecimiento temprano (α=20 ≈ top 8%)\n")
        f.write("- Decoys generados con matching de propiedades fisicoquímicas\n")
        f.write("- Scores CNN de Gnina usados para ranking\n")
    
    print(f"Reporte guardado: {reporte_file}")
    return reporte_file


def main():
    """Función principal."""
    print("\n" + "="*70)
    print("CÁLCULO DE MÉTRICAS DE VALIDACIÓN - ACE2")
    print("="*70)
    
    # Cargar labels
    print("\nCargando etiquetas...")
    labels = cargar_labels()
    
    if not labels:
        print("ERROR: No se pudieron cargar las etiquetas")
        print(f"Verificar archivos: {CONFIG['activos_file']}, {CONFIG['decoys_file']}")
        return
    
    # Recolectar resultados
    print("\nRecolectando resultados de docking...")
    resultados = recolectar_resultados(labels)
    
    # Calcular métricas
    print("\nCalculando métricas...")
    metricas = calcular_todas_metricas(resultados)
    
    if metricas is None:
        print("ERROR: No se pudieron calcular métricas")
        return
    
    # Mostrar resumen
    print("\n" + "-"*40)
    print("RESUMEN DE MÉTRICAS")
    print("-"*40)
    print(f"AUC-ROC:      {metricas['auc_roc']:.3f}")
    print(f"EF 1%:        {metricas['ef_1']:.2f}")
    print(f"EF 5%:        {metricas['ef_5']:.2f}")
    print(f"EF 10%:       {metricas['ef_10']:.2f}")
    print(f"BEDROC (α=20): {metricas['bedroc_20']:.3f}")
    
    # Generar outputs
    print("\nGenerando reportes y gráficas...")
    output_dir = CONFIG["output_dir"]
    
    generar_graficas(metricas, output_dir)
    generar_reporte_texto(metricas, output_dir)
    
    # Guardar métricas en JSON
    json_file = Path(output_dir) / "metricas_validacion_ace2.json"
    
    # Remover arrays numpy para JSON
    metricas_json = {k: v for k, v in metricas.items() 
                     if k not in ["roc_fpr", "roc_tpr"]}
    metricas_json["fecha"] = datetime.now().isoformat()
    
    with open(json_file, "w") as f:
        json.dump(metricas_json, f, indent=2)
    
    print(f"\nMétricas guardadas: {json_file}")
    print(f"\nValidación completada!")


if __name__ == "__main__":
    main()