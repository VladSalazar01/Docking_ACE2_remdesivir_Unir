#!/usr/bin/env python3
"""
Análisis de Resultados DOCK6 Benchmark
Calcula RMSD para cada run y genera estadísticas

TFM: Cribado virtual basado en la estructura de la proteína ACE2
"""

import os
import math
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN
# =============================================================================
BASE_DIR = "/mnt/g/Dock6_benchmark_re"
REFERENCE_FILE = "ligand_ref.pdb"
NUM_RUNS = 5

# =============================================================================
# FUNCIONES
# =============================================================================

def read_pdb_coords(filename):
    """Lee coordenadas de átomos pesados desde archivo PDB"""
    coords = []
    atom_names = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                atom_name = line[12:16].strip()
                # Solo átomos pesados (excluir hidrógenos)
                if not atom_name.startswith('H') and not atom_name.startswith('1H') \
                   and not atom_name.startswith('2H') and not atom_name.startswith('3H'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                    atom_names.append(atom_name)
    
    return coords, atom_names


def read_mol2_coords(filename):
    """Lee coordenadas de átomos pesados desde archivo MOL2"""
    coords = []
    atom_names = []
    in_atoms = False
    
    with open(filename, 'r') as f:
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atoms = True
                continue
            elif '@<TRIPOS>' in line:
                in_atoms = False
                continue
            
            if in_atoms:
                parts = line.split()
                if len(parts) >= 6:
                    try:
                        atom_name = parts[1]
                        x = float(parts[2])
                        y = float(parts[3])
                        z = float(parts[4])
                        atom_type = parts[5]
                        
                        # Solo átomos pesados
                        if not atom_type.startswith('H'):
                            coords.append((x, y, z))
                            atom_names.append(atom_name)
                    except ValueError:
                        pass
    
    return coords, atom_names


def calculate_rmsd(coords1, coords2):
    """Calcula RMSD entre dos conjuntos de coordenadas"""
    if len(coords1) != len(coords2):
        return None
    
    n = len(coords1)
    sum_sq = 0.0
    
    for (x1, y1, z1), (x2, y2, z2) in zip(coords1, coords2):
        sum_sq += (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
    
    rmsd = math.sqrt(sum_sq / n)
    return rmsd


def parse_dock_output(filename):
    """Extrae información del archivo de salida de DOCK6"""
    info = {
        'score': None,
        'vdw_energy': None,
        'es_energy': None,
        'time': None,
        'orientations': None,
        'conformations': None
    }
    
    if not os.path.exists(filename):
        return info
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # Buscar valores
    import re
    
    score_match = re.search(r'Continuous_Score:\s+([-\d.]+)', content)
    if score_match:
        info['score'] = float(score_match.group(1))
    
    vdw_match = re.search(r'Continuous_vdw_energy:\s+([-\d.]+)', content)
    if vdw_match:
        info['vdw_energy'] = float(vdw_match.group(1))
    
    es_match = re.search(r'Continuous_es_energy:\s+([-\d.]+)', content)
    if es_match:
        info['es_energy'] = float(es_match.group(1))
    
    time_match = re.search(r'Elapsed time for docking:\s+([\d.]+)', content)
    if time_match:
        info['time'] = float(time_match.group(1))
    
    orient_match = re.search(r'Orientations:\s+(\d+)', content)
    if orient_match:
        info['orientations'] = int(orient_match.group(1))
    
    conf_match = re.search(r'Conformations:\s+(\d+)', content)
    if conf_match:
        info['conformations'] = int(conf_match.group(1))
    
    return info


def main():
    """Función principal de análisis"""
    
    print("=" * 70)
    print("ANÁLISIS DE RESULTADOS - DOCK6 BENCHMARK")
    print("=" * 70)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Directorio: {BASE_DIR}")
    print("=" * 70)
    
    os.chdir(BASE_DIR)
    
    # Leer coordenadas de referencia
    print(f"\nLeyendo referencia: {REFERENCE_FILE}")
    ref_coords, ref_atoms = read_pdb_coords(REFERENCE_FILE)
    print(f"  Átomos pesados: {len(ref_coords)}")
    
    # Analizar cada run
    results = []
    
    print(f"\nAnalizando {NUM_RUNS} runs...")
    print("-" * 70)
    print(f"{'Run':<6} {'RMSD (Å)':<12} {'Score':<14} {'VDW':<12} {'ES':<12} {'Tiempo (s)':<12}")
    print("-" * 70)
    
    for i in range(1, NUM_RUNS + 1):
        mol2_file = f"output_run{i}_scored.mol2"
        out_file = f"dock_run{i}.out"
        
        if not os.path.exists(mol2_file):
            print(f"{i:<6} {'N/A':<12} {'Archivo no encontrado'}")
            continue
        
        # Leer pose predicha
        pred_coords, pred_atoms = read_mol2_coords(mol2_file)
        
        # Calcular RMSD
        rmsd = calculate_rmsd(ref_coords, pred_coords)
        
        # Leer info del output
        info = parse_dock_output(out_file)
        
        if rmsd is not None:
            results.append({
                'run': i,
                'rmsd': rmsd,
                'score': info['score'],
                'vdw': info['vdw_energy'],
                'es': info['es_energy'],
                'time': info['time']
            })
            
            score_str = f"{info['score']:.2f}" if info['score'] else "N/A"
            vdw_str = f"{info['vdw_energy']:.2f}" if info['vdw_energy'] else "N/A"
            es_str = f"{info['es_energy']:.2f}" if info['es_energy'] else "N/A"
            time_str = f"{info['time']:.1f}" if info['time'] else "N/A"
            
            print(f"{i:<6} {rmsd:<12.3f} {score_str:<14} {vdw_str:<12} {es_str:<12} {time_str:<12}")
        else:
            print(f"{i:<6} {'ERROR':<12} {'Átomos no coinciden'}")
    
    # Estadísticas
    if results:
        print("\n" + "=" * 70)
        print("ESTADÍSTICAS")
        print("=" * 70)
        
        rmsds = [r['rmsd'] for r in results]
        scores = [r['score'] for r in results if r['score']]
        times = [r['time'] for r in results if r['time']]
        
        # RMSD
        rmsd_avg = sum(rmsds) / len(rmsds)
        rmsd_min = min(rmsds)
        rmsd_max = max(rmsds)
        rmsd_std = math.sqrt(sum((x - rmsd_avg)**2 for x in rmsds) / len(rmsds))
        
        print(f"\nRMSD (Å):")
        print(f"  Promedio:          {rmsd_avg:.3f}")
        print(f"  Desv. Estándar:    {rmsd_std:.3f}")
        print(f"  Mínimo:            {rmsd_min:.3f}")
        print(f"  Máximo:            {rmsd_max:.3f}")
        print(f"  Rango:             {rmsd_min:.3f} - {rmsd_max:.3f}")
        
        # Score
        if scores:
            score_avg = sum(scores) / len(scores)
            print(f"\nContinuous Score (kcal/mol):")
            print(f"  Promedio:          {score_avg:.2f}")
        
        # Tiempo
        if times:
            time_avg = sum(times) / len(times)
            print(f"\nTiempo por run:")
            print(f"  Promedio:          {time_avg:.1f} segundos ({time_avg/60:.1f} min)")
        
        # Evaluación
        print("\n" + "=" * 70)
        print("EVALUACIÓN")
        print("=" * 70)
        
        if rmsd_avg < 2.0:
            print(f"✓ RMSD promedio ({rmsd_avg:.3f} Å) < 2.0 Å: EXCELENTE")
        elif rmsd_avg < 3.0:
            print(f"○ RMSD promedio ({rmsd_avg:.3f} Å) < 3.0 Å: BUENO")
        else:
            print(f"✗ RMSD promedio ({rmsd_avg:.3f} Å) >= 3.0 Å: POBRE")
        
        # Guardar resultados en CSV
        csv_file = "dock6_benchmark_results.csv"
        with open(csv_file, 'w') as f:
            f.write("Run,RMSD_A,Continuous_Score,VDW_Energy,ES_Energy,Time_s\n")
            for r in results:
                f.write(f"{r['run']},{r['rmsd']:.3f},{r['score'] or ''},{r['vdw'] or ''},{r['es'] or ''},{r['time'] or ''}\n")
        
        print(f"\nResultados guardados en: {csv_file}")
    
    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()
