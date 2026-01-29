#!/usr/bin/env python3
"""Verificar contenido de logs de Gnina."""
from pathlib import Path

salidas_dir = Path("resultados_validacion_ace2/salidas_docking")

# Buscar archivos .log
logs = list(salidas_dir.glob("*.log"))
print(f"Total archivos .log: {len(logs)}")

if logs:
    # Mostrar contenido de un log
    log = logs[0]
    print(f"\n--- Contenido de {log.name} ---")
    with open(log, "r") as f:
        print(f.read())
