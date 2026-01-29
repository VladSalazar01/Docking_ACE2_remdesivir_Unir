#!/usr/bin/env python3
"""Diagnóstico de scores: activos vs decoys."""
import re
from pathlib import Path

# Cargar labels
activos = set()
decoys = set()

with open("dataset_ace2/validacion/activos.smi") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            activos.add(parts[1])

with open("dataset_ace2/validacion/decoys.smi") as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 2:
            decoys.add(parts[1])

print(f"Activos: {len(activos)}")
print(f"Decoys: {len(decoys)}")

# Extraer scores de logs
salidas_dir = Path("resultados_validacion_ace2/salidas_docking")

scores_activos = []
scores_decoys = []

for log in salidas_dir.glob("*.log"):
    nombre = log.stem
    
    with open(log, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    
    # Extraer CNN score
    pattern = r'^\s*1\s+(-?[\d.]+)\s+(-?[\d.]+)\s+([\d.]+)\s+(-?[\d.]+)'
    match = re.search(pattern, content, re.MULTILINE)
    
    if match:
        cnn_score = float(match.group(3))
        cnn_affinity = float(match.group(4))
        vina = float(match.group(1))
        
        if nombre in activos:
            scores_activos.append((nombre, cnn_score, cnn_affinity, vina))
        elif nombre in decoys:
            scores_decoys.append((nombre, cnn_score, cnn_affinity, vina))

print("\n" + "="*60)
print("ACTIVOS (ordenados por CNN score)")
print("="*60)
print(f"{'Nombre':<20} {'CNN_score':>10} {'CNN_aff':>10} {'Vina':>10}")
print("-"*60)
for nombre, cnn, aff, vina in sorted(scores_activos, key=lambda x: -x[1]):
    print(f"{nombre:<20} {cnn:>10.4f} {aff:>10.2f} {vina:>10.2f}")

print(f"\nMedia CNN score activos: {sum(s[1] for s in scores_activos)/len(scores_activos):.4f}")

print("\n" + "="*60)
print("DECOYS (top 10 por CNN score)")
print("="*60)
print(f"{'Nombre':<20} {'CNN_score':>10} {'CNN_aff':>10} {'Vina':>10}")
print("-"*60)
for nombre, cnn, aff, vina in sorted(scores_decoys, key=lambda x: -x[1])[:10]:
    print(f"{nombre:<20} {cnn:>10.4f} {aff:>10.2f} {vina:>10.2f}")

print(f"\nMedia CNN score decoys: {sum(s[1] for s in scores_decoys)/len(scores_decoys):.4f}")

# Estadísticas
print("\n" + "="*60)
print("RESUMEN ESTADÍSTICO")
print("="*60)
media_act = sum(s[1] for s in scores_activos)/len(scores_activos)
media_dec = sum(s[1] for s in scores_decoys)/len(scores_decoys)
print(f"Media CNN score activos: {media_act:.4f}")
print(f"Media CNN score decoys:  {media_dec:.4f}")
print(f"Diferencia: {media_act - media_dec:.4f}")

if media_act < media_dec:
    print("\n⚠️  PROBLEMA: Los decoys tienen scores más altos que los activos!")
    print("   Esto explica el AUC-ROC < 0.5")
