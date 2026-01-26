#!/usr/bin/env python3
"""
CÁLCULO DE MÉTRICAS DE VALIDACIÓN RETROSPECTIVA
Carpeta de trabajo: G:\Cribado_VC_ACE_2\validacion_ACE1
"""

import csv
import sys
import math
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(r"G:\Cribado_VC_ACE_2\validacion_ACE1")
RESULTS_FILE = WORK_DIR / "results" / "docked_results.sdf"
LABELS_FILE = WORK_DIR / "subset" / "labels.csv"

def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}")

def parse_gnina_results(sdf_file):
    results = {}
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
    i = 0
    current_name = None
    while i < len(lines):
        line = lines[i].strip()
        if i == 0 or (i > 0 and lines[i-1].strip() == '$$$$'):
            if line and not line.startswith('>') and not line.startswith('M  '):
                current_name = line
        if line == '> <CNNaffinity>':
            i += 1
            if i < len(lines):
                try:
                    score = float(lines[i].strip())
                    if current_name:
                        results[current_name] = score
                except ValueError:
                    pass
        if line == '$$$$':
            current_name = None
        i += 1
    return results

def load_labels(csv_file):
    labels = {}
    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            labels[row['mol_name']] = int(row['is_active'])
    return labels

def calculate_auc_roc(scores_labels):
    sorted_data = sorted(scores_labels, key=lambda x: x[0], reverse=True)
    n_pos = sum(1 for _, label in sorted_data if label == 1)
    n_neg = len(sorted_data) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    auc = 0.0
    for i, (score_i, label_i) in enumerate(sorted_data):
        if label_i == 0:
            auc += sum(1 for _, label_j in sorted_data[:i] if label_j == 1)
    auc /= (n_pos * n_neg)
    return auc

def calculate_enrichment(scores_labels, fraction):
    sorted_data = sorted(scores_labels, key=lambda x: x[0], reverse=True)
    n_total = len(sorted_data)
    n_actives_total = sum(1 for _, label in sorted_data if label == 1)
    if n_actives_total == 0:
        return 0.0
    n_top = max(1, int(n_total * fraction))
    n_actives_top = sum(1 for _, label in sorted_data[:n_top] if label == 1)
    expected = n_actives_total * fraction
    if expected == 0:
        return 0.0
    ef = n_actives_top / expected
    return ef

def calculate_bedroc(scores_labels, alpha=80.5):
    sorted_data = sorted(scores_labels, key=lambda x: x[0], reverse=True)
    n = len(sorted_data)
    n_actives = sum(1 for _, label in sorted_data if label == 1)
    if n_actives == 0 or n_actives == n:
        return 0.0
    sum_exp = 0.0
    for i, (score, label) in enumerate(sorted_data):
        if label == 1:
            rank = (i + 1) / n
            sum_exp += math.exp(-alpha * rank)
    ra = n_actives / n
    random_sum = ra * n * (1 - math.exp(-alpha)) / (1 - math.exp(-alpha / n)) / n
    ri_max = (1 - math.exp(-alpha * ra)) / (1 - math.exp(-alpha / n))
    ri_min = ra * (1 - math.exp(-alpha)) / (1 - math.exp(-alpha / n))
    if ri_max == ri_min:
        return 0.0
    bedroc = (sum_exp / n_actives - ri_min / n_actives) / (ri_max / n_actives - ri_min / n_actives)
    bedroc = max(0, min(1, bedroc))
    return bedroc

def main():
    log("="*60)
    log("METRICAS DE VALIDACION RETROSPECTIVA")
    log("TFM: Cribado Virtual ACE2")
    log("="*60)
    if not RESULTS_FILE.exists():
        log(f"ERROR: Resultados no encontrados: {RESULTS_FILE}")
        sys.exit(1)
    if not LABELS_FILE.exists():
        log(f"ERROR: Etiquetas no encontradas: {LABELS_FILE}")
        sys.exit(1)
    log("")
    log(f"Cargando resultados: {RESULTS_FILE.name}")
    scores = parse_gnina_results(RESULTS_FILE)
    log(f"  {len(scores)} moleculas con scores")
    log(f"Cargando etiquetas: {LABELS_FILE.name}")
    labels = load_labels(LABELS_FILE)
    log(f"  {len(labels)} etiquetas")
    scores_labels = []
    missing = 0
    for name, label in labels.items():
        if name in scores:
            scores_labels.append((scores[name], label))
        else:
            missing += 1
    if missing > 0:
        log(f"  ADVERTENCIA: {missing} moleculas sin score")
    n_actives = sum(1 for _, l in scores_labels if l == 1)
    n_decoys = len(scores_labels) - n_actives
    log("")
    log(f"Datos combinados: {len(scores_labels)} moleculas")
    log(f"  Activos: {n_actives}")
    log(f"  Decoys: {n_decoys}")
    log("")
    log("="*60)
    log("RESULTADOS")
    log("="*60)
    auc = calculate_auc_roc(scores_labels)
    ef1 = calculate_enrichment(scores_labels, 0.01)
    ef5 = calculate_enrichment(scores_labels, 0.05)
    bedroc = calculate_bedroc(scores_labels, alpha=80.5)
    log("")
    log(f"AUC-ROC:          {auc:.4f}    {'PASS' if auc >= 0.70 else 'FAIL'} (umbral >= 0.70)")
    log(f"EF 1%:            {ef1:.2f}     {'PASS' if ef1 >= 10 else 'FAIL'} (umbral >= 10)")
    log(f"EF 5%:            {ef5:.2f}     {'PASS' if ef5 >= 5 else 'FAIL'} (umbral >= 5)")
    log(f"BEDROC (a=80.5):  {bedroc:.4f}    {'PASS' if bedroc >= 0.30 else 'FAIL'} (umbral >= 0.30)")
    log("")
    log("="*60)
    log("EVALUACION GLOBAL")
    log("="*60)
    if auc >= 0.85 and ef1 >= 20:
        log(">>> OPTIMO: Proceder con screening prospectivo")
        evaluacion = "OPTIMO"
    elif auc >= 0.70 and ef1 >= 10:
        log(">>> ACEPTABLE: Proceder con precaucion")
        evaluacion = "ACEPTABLE"
    else:
        log(">>> SUBOPTIMO: Requiere optimizacion del protocolo")
        evaluacion = "SUBOPTIMO"
    log("")
    log("-"*60)
    log("TOP 10 PREDICCIONES")
    log("-"*60)
    name_score_label = []
    for name, label in labels.items():
        if name in scores:
            name_score_label.append((name, scores[name], label))
    name_score_label.sort(key=lambda x: x[1], reverse=True)
    log(f"{'Rank':<6} {'Nombre':<20} {'Score':<12} {'Real':<10}")
    log("-"*60)
    for i, (name, score, label) in enumerate(name_score_label[:10], 1):
        label_str = "ACTIVO" if label == 1 else "decoy"
        log(f"{i:<6} {name[:20]:<20} {score:<12.4f} {label_str:<10}")
    output_file = WORK_DIR / "results" / "metricas_validacion.txt"
    with open(output_file, 'w') as f:
        f.write("METRICAS DE VALIDACION RETROSPECTIVA\n")
        f.write("="*50 + "\n\n")
        f.write(f"Activos: {n_actives}\n")
        f.write(f"Decoys: {n_decoys}\n")
        f.write(f"Total: {len(scores_labels)}\n\n")
        f.write(f"AUC-ROC: {auc:.4f}\n")
        f.write(f"EF 1%: {ef1:.2f}\n")
        f.write(f"EF 5%: {ef5:.2f}\n")
        f.write(f"BEDROC (alpha=80.5): {bedroc:.4f}\n\n")
        f.write(f"Evaluacion: {evaluacion}\n")
    log("")
    log(f"Resultados guardados en: {output_file}")

if __name__ == "__main__":
    main()