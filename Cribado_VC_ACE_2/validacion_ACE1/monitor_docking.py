#!/usr/bin/env python3
"""
Monitor de progreso de docking - Métricas en tiempo real
Carpeta: G:\Cribado_VC_ACE_2\validacion_ACE1
"""

import time
import csv
from pathlib import Path
from datetime import datetime

WORK_DIR = Path(r"G:\Cribado_VC_ACE_2\validacion_ACE1")
RESULTS_FILE = WORK_DIR / "results" / "docked_results.sdf"
LABELS_FILE = WORK_DIR / "subset" / "labels.csv"
INTERVALO = 60

def load_labels():
    labels = {}
    with open(LABELS_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            labels[row['mol_name']] = int(row['is_active'])
    return labels

def parse_partial_results(sdf_file):
    results = {}
    try:
        with open(sdf_file, 'r') as f:
            content = f.read()
        blocks = content.split('$$$$')
        for block in blocks:
            if not block.strip():
                continue
            lines = block.strip().split('\n')
            name = lines[0].strip() if lines else None
            score = None
            for i, line in enumerate(lines):
                if '<CNNaffinity>' in line and i+1 < len(lines):
                    try:
                        score = float(lines[i+1].strip())
                    except:
                        pass
            if name and score:
                results[name] = score
    except:
        pass
    return results

def calc_auc(scores_labels):
    if len(scores_labels) < 10:
        return None
    sorted_data = sorted(scores_labels, key=lambda x: x[0], reverse=True)
    n_pos = sum(1 for _, l in sorted_data if l == 1)
    n_neg = len(sorted_data) - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.5
    auc = sum(1 for i, (_, l1) in enumerate(sorted_data) if l1 == 0 
              for _, l2 in sorted_data[:i] if l2 == 1)
    return auc / (n_pos * n_neg)

def calc_ef1(scores_labels):
    if len(scores_labels) < 100:
        return None
    sorted_data = sorted(scores_labels, key=lambda x: x[0], reverse=True)
    n_total = len(sorted_data)
    n_actives = sum(1 for _, l in sorted_data if l == 1)
    if n_actives == 0:
        return 0
    n_top = max(1, int(n_total * 0.01))
    actives_top = sum(1 for _, l in sorted_data[:n_top] if l == 1)
    return actives_top / (n_actives * 0.01)

def main():
    print("="*60)
    print("MONITOR DE DOCKING EN TIEMPO REAL")
    print("="*60)
    print(f"Archivo: {RESULTS_FILE}")
    print(f"Intervalo: {INTERVALO}s")
    print("Ctrl+C para detener")
    print("="*60)
    labels = load_labels()
    n_total = len(labels)
    n_activos_total = sum(labels.values())
    print(f"Esperados: {n_total} mol ({n_activos_total} activos)")
    print("-"*60)
    historial = []
    while True:
        try:
            scores = parse_partial_results(RESULTS_FILE)
            n_proc = len(scores)
            if n_proc == 0:
                print(f"[{datetime.now().strftime('%H:%M:%S')}] Esperando resultados...")
                time.sleep(INTERVALO)
                continue
            scores_labels = [(scores[n], labels[n]) for n in labels if n in scores]
            n_act = sum(1 for _, l in scores_labels if l == 1)
            auc = calc_auc(scores_labels)
            ef1 = calc_ef1(scores_labels)
            pct = n_proc / n_total * 100
            ts = datetime.now().strftime('%H:%M:%S')
            auc_str = f"{auc:.3f}" if auc else "---"
            ef1_str = f"{ef1:.1f}" if ef1 else "---"
            tendencia = ""
            if auc and len(historial) > 0:
                diff = auc - historial[-1]
                tendencia = "↑" if diff > 0.01 else "↓" if diff < -0.01 else "→"
            if auc:
                historial.append(auc)
            status = "OK" if auc and auc >= 0.70 else "BAJO" if auc else ""
            print(f"[{ts}] {n_proc:4d}/{n_total} ({pct:5.1f}%) | "
                  f"Act:{n_act:2d} | AUC:{auc_str} {tendencia} | EF1%:{ef1_str} | {status}")
            if n_proc >= n_total:
                print("-"*60)
                print("DOCKING COMPLETADO")
                break
            time.sleep(INTERVALO)
        except KeyboardInterrupt:
            print("\nMonitor detenido")
            break

if __name__ == "__main__":
    main()