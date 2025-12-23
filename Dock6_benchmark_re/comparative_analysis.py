#!/usr/bin/env python3
"""
Análisis Comparativo: DOCK6 vs Gnina
Benchmark de Re-docking - Estructura 9FMM (Fragment 1, Ligando A1IDX)

TFM: Cribado virtual basado en la estructura de la proteína ACE2
Autor: [Tu nombre]
Fecha: Diciembre 2025
"""

import os
from datetime import datetime

# =============================================================================
# DATOS DEL BENCHMARK
# =============================================================================

# Configuración común del sistema de prueba
SYSTEM_CONFIG = {
    'pdb_id': '9FMM',
    'ligand_id': 'A1IDX',
    'fragment': 'Fragment 1 (Cadena A, residuo 801)',
    'binding_site_center': (16.48, 15.24, 25.54),
    'box_size': (20, 20, 20),  # Angstroms
    'heavy_atoms_ligand': 28,
    'receptor_atoms_full': 11022,
}

# Configuración del hardware
HARDWARE_CONFIG = {
    'cpu': 'AMD Ryzen Threadripper 3690x (24 cores / 48 threads)',
    'ram': '128 GB DDR4',
    'gpu': '2x NVIDIA RTX 3060 (12 GB VRAM cada una)',
    'os': 'Windows 10 x64 + WSL2 Ubuntu 22.04',
    'storage': 'SSD NVMe',
}

# Resultados de Gnina (Fragment 1)
GNINA_CONFIG = {
    'version': '1.3.1',
    'scoring': 'CNN rescore',
    'exhaustiveness': 32,
    'num_modes': 9,
    'cpu_threads': 48,
    'gpus': '0,1 (Dual RTX 3060)',
    'receptor_format': 'PDB',
    'ligand_format': 'PDB',
    'receptor_atoms': 11022,  # Receptor completo
    'num_runs': 17,
}

GNINA_RESULTS = {
    'rmsd_values': [1.496, 1.495, 1.499, 1.498, 1.497, 1.496, 1.495, 1.498, 
                   1.497, 1.496, 1.495, 1.499, 1.498, 1.497, 1.496, 1.495, 1.498],
    'time_per_run_seconds': 27,  # Aproximado con GPU
    'best_affinity_kcal': -9.43,
    'cnn_pose_score': 0.9706,
}

# Resultados de DOCK6
DOCK6_CONFIG = {
    'version': '6.13',
    'scoring': 'Continuous Score',
    'conformer_search': 'Flexible (anchor-and-grow)',
    'max_orientations': 1000,
    'pruning_max_orients': 1000,
    'simplex_max_iterations': 1000,
    'receptor_format': 'MOL2',
    'ligand_format': 'MOL2',
    'receptor_atoms': 1255,  # Receptor recortado (15 Å del ligando)
    'num_spheres': 59,
    'num_runs': 5,
    'parallelization': 'Single-thread (secuencial)',
}

DOCK6_RESULTS = {
    'rmsd_values': [0.485, 0.467, 0.468, 0.483, 0.462],
    'scores': [-37.81, -37.79, -37.81, -37.80, -37.79],
    'vdw_energies': [-38.11, -38.09, -38.10, -38.08, -38.09],
    'es_energies': [0.30, 0.30, 0.29, 0.29, 0.30],
    'times_seconds': [1286.3, 1313.1, 1357.3, 1337.4, 1523.6],
}

# =============================================================================
# FUNCIONES DE ANÁLISIS
# =============================================================================

def calculate_stats(values):
    """Calcula estadísticas básicas de una lista de valores"""
    n = len(values)
    mean = sum(values) / n
    variance = sum((x - mean) ** 2 for x in values) / n
    std = variance ** 0.5
    return {
        'n': n,
        'mean': mean,
        'std': std,
        'min': min(values),
        'max': max(values),
        'range': max(values) - min(values),
    }


def generate_report():
    """Genera el reporte comparativo completo"""
    
    # Calcular estadísticas
    gnina_stats = calculate_stats(GNINA_RESULTS['rmsd_values'])
    dock6_stats = calculate_stats(DOCK6_RESULTS['rmsd_values'])
    dock6_time_stats = calculate_stats(DOCK6_RESULTS['times_seconds'])
    dock6_score_stats = calculate_stats(DOCK6_RESULTS['scores'])
    
    report = []
    
    # Encabezado
    report.append("=" * 80)
    report.append("ANÁLISIS COMPARATIVO: DOCK6 vs GNINA")
    report.append("Benchmark de Re-docking - Estructura 9FMM")
    report.append("=" * 80)
    report.append("Fecha de análisis: {}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    report.append("")
    
    # Sistema de prueba
    report.append("-" * 80)
    report.append("1. SISTEMA DE PRUEBA")
    report.append("-" * 80)
    report.append("  Estructura PDB:        {}".format(SYSTEM_CONFIG['pdb_id']))
    report.append("  Ligando:               {}".format(SYSTEM_CONFIG['ligand_id']))
    report.append("  Sitio de unión:        {}".format(SYSTEM_CONFIG['fragment']))
    center = SYSTEM_CONFIG['binding_site_center']
    report.append("  Centro del sitio:      ({}, {}, {}) Å".format(center[0], center[1], center[2]))
    box = SYSTEM_CONFIG['box_size']
    report.append("  Tamaño del box:        {} × {} × {} Å".format(box[0], box[1], box[2]))
    report.append("  Átomos pesados ligando: {}".format(SYSTEM_CONFIG['heavy_atoms_ligand']))
    report.append("")
    
    # Hardware
    report.append("-" * 80)
    report.append("2. CONFIGURACIÓN DE HARDWARE")
    report.append("-" * 80)
    report.append("  CPU:                   {}".format(HARDWARE_CONFIG['cpu']))
    report.append("  RAM:                   {}".format(HARDWARE_CONFIG['ram']))
    report.append("  GPU:                   {}".format(HARDWARE_CONFIG['gpu']))
    report.append("  Sistema Operativo:     {}".format(HARDWARE_CONFIG['os']))
    report.append("")
    
    # Configuración de programas
    report.append("-" * 80)
    report.append("3. CONFIGURACIÓN DE PROGRAMAS")
    report.append("-" * 80)
    report.append("")
    report.append("  GNINA:")
    report.append("    Versión:             {}".format(GNINA_CONFIG['version']))
    report.append("    Método de scoring:   {}".format(GNINA_CONFIG['scoring']))
    report.append("    Exhaustiveness:      {}".format(GNINA_CONFIG['exhaustiveness']))
    report.append("    Num modes:           {}".format(GNINA_CONFIG['num_modes']))
    report.append("    Threads CPU:         {}".format(GNINA_CONFIG['cpu_threads']))
    report.append("    GPUs:                {}".format(GNINA_CONFIG['gpus']))
    report.append("    Formato receptor:    {}".format(GNINA_CONFIG['receptor_format']))
    report.append("    Átomos receptor:     {} (completo)".format(GNINA_CONFIG['receptor_atoms']))
    report.append("    Número de runs:      {}".format(GNINA_CONFIG['num_runs']))
    report.append("")
    report.append("  DOCK6:")
    report.append("    Versión:             {}".format(DOCK6_CONFIG['version']))
    report.append("    Método de scoring:   {}".format(DOCK6_CONFIG['scoring']))
    report.append("    Búsqueda conformac.: {}".format(DOCK6_CONFIG['conformer_search']))
    report.append("    Max orientaciones:   {}".format(DOCK6_CONFIG['max_orientations']))
    report.append("    Formato receptor:    {}".format(DOCK6_CONFIG['receptor_format']))
    report.append("    Átomos receptor:     {} (recortado 15 Å)".format(DOCK6_CONFIG['receptor_atoms']))
    report.append("    Número de esferas:   {}".format(DOCK6_CONFIG['num_spheres']))
    report.append("    Número de runs:      {}".format(DOCK6_CONFIG['num_runs']))
    report.append("    Paralelización:      {}".format(DOCK6_CONFIG['parallelization']))
    report.append("")
    
    # Resultados
    report.append("-" * 80)
    report.append("4. RESULTADOS DE RE-DOCKING")
    report.append("-" * 80)
    report.append("")
    report.append("  4.1 RMSD (Root Mean Square Deviation)")
    report.append("")
    report.append("  {:<25} {:<20} {:<20}".format("Métrica", "GNINA", "DOCK6"))
    report.append("  {} {} {}".format("-"*25, "-"*20, "-"*20))
    report.append("  {:<25} {:<20} {:<20}".format("Número de runs", gnina_stats['n'], dock6_stats['n']))
    report.append("  {:<25} {:<20.3f} {:<20.3f}".format("RMSD promedio (Å)", gnina_stats['mean'], dock6_stats['mean']))
    report.append("  {:<25} {:<20.3f} {:<20.3f}".format("Desviación estándar (Å)", gnina_stats['std'], dock6_stats['std']))
    report.append("  {:<25} {:<20.3f} {:<20.3f}".format("RMSD mínimo (Å)", gnina_stats['min'], dock6_stats['min']))
    report.append("  {:<25} {:<20.3f} {:<20.3f}".format("RMSD máximo (Å)", gnina_stats['max'], dock6_stats['max']))
    report.append("  {:<25} {:<20.3f} {:<20.3f}".format("Rango (Å)", gnina_stats['range'], dock6_stats['range']))
    report.append("")
    
    # Mejora porcentual
    improvement = ((gnina_stats['mean'] - dock6_stats['mean']) / gnina_stats['mean']) * 100
    report.append("  Mejora de DOCK6 vs GNINA: {:.1f}% en RMSD promedio".format(improvement))
    report.append("")
    
    # Tiempo computacional
    report.append("  4.2 TIEMPO COMPUTACIONAL")
    report.append("")
    gnina_total = GNINA_RESULTS['time_per_run_seconds'] * GNINA_CONFIG['num_runs']
    dock6_total = sum(DOCK6_RESULTS['times_seconds'])
    
    dock6_time_min = dock6_time_stats['mean'] / 60
    
    report.append("  {:<30} {:<20} {:<20}".format("Métrica", "GNINA", "DOCK6"))
    report.append("  {} {} {}".format("-"*30, "-"*20, "-"*20))
    report.append("  {:<30} {:<20} {:<20}".format("Tiempo por run", "~27 seg (GPU)", "{:.1f} min (CPU)".format(dock6_time_min)))
    report.append("  {:<30} {:<20} {:<20}".format("Tiempo total benchmark", "~{:.1f} min".format(gnina_total/60), "{:.1f} min".format(dock6_total/60)))
    report.append("  {:<30} {:<20} {:<20}".format("Aceleración por hardware", "GPU (2x RTX 3060)", "CPU single-thread"))
    report.append("")
    
    # Scoring
    report.append("  4.3 SCORES DE AFINIDAD")
    report.append("")
    report.append("  {:<30} {:<20} {:<20}".format("Métrica", "GNINA", "DOCK6"))
    report.append("  {} {} {}".format("-"*30, "-"*20, "-"*20))
    gnina_affinity = "{:.2f} kcal/mol".format(GNINA_RESULTS['best_affinity_kcal'])
    dock6_affinity = "{:.2f} kcal/mol".format(dock6_score_stats['mean'])
    report.append("  {:<30} {:<20} {:<20}".format("Mejor afinidad", gnina_affinity, dock6_affinity))
    cnn_score = "{:.4f}".format(GNINA_RESULTS['cnn_pose_score'])
    report.append("  {:<30} {:<20} {:<20}".format("CNN pose score", cnn_score, "N/A"))
    report.append("")
    
    # Evaluación según criterios estándar
    report.append("-" * 80)
    report.append("5. EVALUACIÓN SEGÚN CRITERIOS ESTÁNDAR")
    report.append("-" * 80)
    report.append("")
    report.append("  Criterios de interpretación de RMSD (literatura):")
    report.append("    - RMSD < 2.0 Å:  Excelente reproducción de la pose experimental")
    report.append("    - RMSD 2.0-3.0 Å: Buena predicción")
    report.append("    - RMSD 3.0-5.0 Å: Aceptable, con desviaciones significativas")
    report.append("    - RMSD > 5.0 Å:  Predicción pobre, pose incorrecta")
    report.append("")
    
    # Evaluación GNINA
    if gnina_stats['mean'] < 2.0:
        gnina_eval = "EXCELENTE"
    elif gnina_stats['mean'] < 3.0:
        gnina_eval = "BUENO"
    elif gnina_stats['mean'] < 5.0:
        gnina_eval = "ACEPTABLE"
    else:
        gnina_eval = "POBRE"
    
    # Evaluación DOCK6
    if dock6_stats['mean'] < 2.0:
        dock6_eval = "EXCELENTE"
    elif dock6_stats['mean'] < 3.0:
        dock6_eval = "BUENO"
    elif dock6_stats['mean'] < 5.0:
        dock6_eval = "ACEPTABLE"
    else:
        dock6_eval = "POBRE"
    
    report.append("  Evaluación GNINA:  {} (RMSD = {:.3f} Å)".format(gnina_eval, gnina_stats['mean']))
    report.append("  Evaluación DOCK6:  {} (RMSD = {:.3f} Å)".format(dock6_eval, dock6_stats['mean']))
    report.append("")
    
    # Diferencias metodológicas
    report.append("-" * 80)
    report.append("6. DIFERENCIAS METODOLÓGICAS CLAVE")
    report.append("-" * 80)
    report.append("")
    report.append("  Aspecto                    GNINA                      DOCK6")
    report.append("  " + "-"*76)
    report.append("  Algoritmo de búsqueda      Monte Carlo + Local Opt    Anchor-and-Grow")
    report.append("  Función de scoring         CNN + Vina empírico        Continuous (VDW + ES)")
    report.append("  Definición del sitio       Box (coordenadas XYZ)      Esferas (superficie)")
    report.append("  Aceleración GPU            Sí (CUDA)                  No")
    report.append("  Preparación receptor       PDB directo                MOL2 con cargas")
    report.append("  Flexibilidad ligando       Torsiones predefinidas     Fragmentación dinámica")
    report.append("")
    
    # Consideraciones
    report.append("-" * 80)
    report.append("7. CONSIDERACIONES PARA LA INTERPRETACIÓN")
    report.append("-" * 80)
    report.append("")
    report.append("  7.1 Ventajas de los resultados de DOCK6:")
    report.append("      - RMSD significativamente menor (0.473 vs 1.497 Å)")
    report.append("      - Alta reproducibilidad (desv. est. = 0.009 Å)")
    report.append("      - No requiere GPU")
    report.append("")
    report.append("  7.2 Ventajas de GNINA:")
    report.append("      - Tiempo de cálculo ~50x más rápido con GPU")
    report.append("      - No requiere recorte del receptor")
    report.append("      - Scoring basado en CNN puede generalizar mejor")
    report.append("")
    report.append("  7.3 Factores que pueden influir en la diferencia de RMSD:")
    report.append("      - DOCK6 usó receptor recortado (1,255 átomos vs 11,022)")
    report.append("      - Diferentes funciones de scoring")
    report.append("      - Diferentes algoritmos de búsqueda conformacional")
    report.append("      - Número de orientaciones/exhaustiveness diferentes")
    report.append("")
    
    # Conclusiones
    report.append("-" * 80)
    report.append("8. CONCLUSIONES")
    report.append("-" * 80)
    report.append("")
    report.append("  1. Ambos programas logran reproducir la pose cristalográfica con")
    report.append("     precisión excelente (RMSD < 2.0 Å), validando su utilidad para")
    report.append("     re-docking en el sistema 9FMM/A1IDX.")
    report.append("")
    report.append("  2. DOCK6 demostró mayor precisión geométrica (RMSD 68% menor),")
    report.append("     aunque a costa de mayor tiempo computacional.")
    report.append("")
    report.append("  3. GNINA ofrece un balance óptimo velocidad/precisión gracias a")
    report.append("     la aceleración GPU y scoring CNN.")
    report.append("")
    report.append("  4. La elección entre programas dependerá del objetivo:")
    report.append("     - Virtual screening de grandes bibliotecas: GNINA (velocidad)")
    report.append("     - Refinamiento de poses específicas: DOCK6 (precisión)")
    report.append("")
    report.append("=" * 80)
    
    return "\n".join(report)


def generate_table_for_tfm():
    """Genera tabla formateada para el documento TFM"""
    
    gnina_stats = calculate_stats(GNINA_RESULTS['rmsd_values'])
    dock6_stats = calculate_stats(DOCK6_RESULTS['rmsd_values'])
    dock6_time_stats = calculate_stats(DOCK6_RESULTS['times_seconds'])
    
    table = []
    table.append("")
    table.append("TABLA: Comparación de resultados de re-docking GNINA vs DOCK6")
    table.append("")
    table.append("| Parámetro | GNINA | DOCK6 |")
    table.append("|-----------|-------|-------|")
    table.append("| Versión | {} | {} |".format(GNINA_CONFIG['version'], DOCK6_CONFIG['version']))
    table.append("| Método de scoring | {} | {} |".format(GNINA_CONFIG['scoring'], DOCK6_CONFIG['scoring']))
    table.append("| Número de runs | {} | {} |".format(GNINA_CONFIG['num_runs'], DOCK6_CONFIG['num_runs']))
    table.append("| RMSD promedio (Å) | {:.3f} | {:.3f} |".format(gnina_stats['mean'], dock6_stats['mean']))
    table.append("| Desv. estándar (Å) | {:.3f} | {:.3f} |".format(gnina_stats['std'], dock6_stats['std']))
    table.append("| Rango RMSD (Å) | {:.3f} - {:.3f} | {:.3f} - {:.3f} |".format(
        gnina_stats['min'], gnina_stats['max'], dock6_stats['min'], dock6_stats['max']))
    table.append("| Tiempo/run | ~27 seg | ~{:.1f} min |".format(dock6_time_stats['mean']/60))
    table.append("| Aceleración | GPU | CPU |")
    table.append("| Evaluación | EXCELENTE | EXCELENTE |")
    table.append("")
    
    return "\n".join(table)


def main():
    """Función principal"""
    
    # Generar reporte completo
    report = generate_report()
    print(report)
    
    # Guardar reporte
    report_file = "comparative_analysis_dock6_gnina.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)
    print("\nReporte guardado en: {}".format(report_file))
    
    # Generar tabla para TFM
    table = generate_table_for_tfm()
    print(table)
    
    table_file = "table_comparison_tfm.md"
    with open(table_file, 'w', encoding='utf-8') as f:
        f.write(table)
    print("Tabla guardada en: {}".format(table_file))


if __name__ == "__main__":
    main()