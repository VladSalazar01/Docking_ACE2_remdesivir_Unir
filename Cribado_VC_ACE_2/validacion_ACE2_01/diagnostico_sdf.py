#!/usr/bin/env python3
"""Diagnóstico de archivos SDF de Gnina."""
from pathlib import Path

salidas_dir = Path("resultados_validacion_ace2/salidas_docking")

# Verificar estructura de un archivo SDF
archivos = list(salidas_dir.glob("*_out.sdf"))
print(f"Total archivos SDF: {len(archivos)}")

con_score = 0
sin_score = []

for sdf in archivos:
    with open(sdf, "r") as f:
        contenido = f.read()
    
    tiene_cnn = "> <CNNscore>" in contenido
    if tiene_cnn:
        con_score += 1
    else:
        sin_score.append(sdf.stem)
        # Mostrar primeras líneas para debug
        if len(sin_score) <= 3:
            print(f"\n--- {sdf.name} (sin CNNscore) ---")
            print(contenido[:500])

print(f"\nCon CNNscore: {con_score}")
print(f"Sin CNNscore: {len(sin_score)}")
if sin_score:
    print(f"Ejemplos sin score: {sin_score[:10]}")

# Mostrar ejemplo de archivo CON score
for sdf in archivos:
    with open(sdf, "r") as f:
        contenido = f.read()
    if "> <CNNscore>" in contenido:
        print(f"\n--- EJEMPLO CON SCORE: {sdf.name} ---")
        # Mostrar solo la parte de propiedades
        idx = contenido.find("> <")
        if idx > 0:
            print(contenido[idx:idx+300])
        break