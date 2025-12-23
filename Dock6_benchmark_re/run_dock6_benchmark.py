#!/usr/bin/env python3
"""
DOCK6 Benchmark - 5 Runs para Re-docking de 9FMM (Fragment 1)
TFM: Cribado virtual basado en la estructura de la proteína ACE2

Autor: [Tu nombre]
Fecha: Diciembre 2025
Sistema: AMD Threadripper 3690x, 128GB RAM, Windows 10 + WSL2 Ubuntu 22.04
"""

import subprocess
import os
import time
from datetime import datetime

# =============================================================================
# CONFIGURACIÓN
# =============================================================================
NUM_RUNS = 5
BASE_DIR = "/mnt/g/Dock6_benchmark_re"
TEMPLATE_FILE = "dock.in"
RECEPTOR_FILE = "receptor_site.mol2"
LIGAND_FILE = "ligand.mol2"
SPHERES_FILE = "selected_spheres.sph"

# =============================================================================
# FUNCIONES
# =============================================================================

def create_dock_input(run_id, seed):
    """Crea archivo de configuración para un run específico"""
    
    input_content = f"""conformer_search_type                                        flex
write_fragment_libraries                                     no
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              5
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               100.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             {LIGAND_FILE}
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               yes
use_rmsd_reference_mol                                       yes
rmsd_reference_filename                                      {LIGAND_FILE}
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           {SPHERES_FILE}
max_orientations                                             1000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      no
multigrid_score_secondary                                    no
dock3.5_score_primary                                        no
dock3.5_score_secondary                                      no
continuous_score_primary                                     yes
cont_score_rec_filename                                      {RECEPTOR_FILE}
cont_score_att_exp                                           6
cont_score_rep_exp                                           12
cont_score_rep_rad_scale                                     1
cont_score_use_dist_dep_dielectric                           yes
cont_score_dielectric                                        4.0
cont_score_vdw_scale                                         1
cont_score_es_scale                                          1
continuous_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          {seed}
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                vdw_AMBER_parm99.defn
flex_defn_file                                               flex.defn
flex_drive_file                                              flex_drive.tbl
ligand_outfile_prefix                                        output_run{run_id}
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
"""
    
    filename = f"dock_run{run_id}.in"
    with open(filename, 'w') as f:
        f.write(input_content)
    
    return filename


def run_dock6(run_id):
    """Ejecuta un run de DOCK6 y retorna tiempo de ejecución"""
    
    input_file = f"dock_run{run_id}.in"
    output_file = f"dock_run{run_id}.out"
    
    print(f"  Run {run_id}: Iniciando...")
    start_time = time.time()
    
    result = subprocess.run(
        ["dock6", "-i", input_file, "-o", output_file],
        capture_output=True,
        text=True
    )
    
    elapsed_time = time.time() - start_time
    
    if result.returncode == 0:
        print(f"  Run {run_id}: Completado en {elapsed_time:.1f} segundos")
    else:
        print(f"  Run {run_id}: ERROR - {result.stderr[:200]}")
    
    return elapsed_time, result.returncode


def main():
    """Función principal del benchmark"""
    
    print("=" * 70)
    print("DOCK6 BENCHMARK - Re-docking 9FMM Fragment 1")
    print("=" * 70)
    print(f"Fecha: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Directorio: {BASE_DIR}")
    print(f"Número de runs: {NUM_RUNS}")
    print("=" * 70)
    
    os.chdir(BASE_DIR)
    
    # Verificar archivos necesarios
    required_files = [RECEPTOR_FILE, LIGAND_FILE, SPHERES_FILE, 
                      "vdw_AMBER_parm99.defn", "flex.defn", "flex_drive.tbl"]
    
    print("\nVerificando archivos...")
    for f in required_files:
        if os.path.exists(f):
            print(f"  ✓ {f}")
        else:
            print(f"  ✗ {f} - NO ENCONTRADO")
            return
    
    # Crear archivos de configuración
    print(f"\nCreando {NUM_RUNS} archivos de configuración...")
    for i in range(1, NUM_RUNS + 1):
        seed = i * 10  # Seeds: 10, 20, 30, 40, 50
        filename = create_dock_input(i, seed)
        print(f"  Creado: {filename} (seed={seed})")
    
    # Ejecutar runs SECUENCIALMENTE (más seguro para memoria)
    print(f"\nEjecutando {NUM_RUNS} runs secuencialmente...")
    print("(Cada run toma ~24 minutos)")
    print("-" * 70)
    
    total_start = time.time()
    results = []
    
    for i in range(1, NUM_RUNS + 1):
        elapsed, returncode = run_dock6(i)
        results.append({
            'run': i,
            'time': elapsed,
            'success': returncode == 0
        })
    
    total_time = time.time() - total_start
    
    # Resumen
    print("\n" + "=" * 70)
    print("RESUMEN DE EJECUCIÓN")
    print("=" * 70)
    
    successful = sum(1 for r in results if r['success'])
    print(f"Runs exitosos: {successful}/{NUM_RUNS}")
    print(f"Tiempo total: {total_time/60:.1f} minutos")
    
    if successful > 0:
        avg_time = sum(r['time'] for r in results if r['success']) / successful
        print(f"Tiempo promedio por run: {avg_time/60:.1f} minutos")
    
    print("\nArchivos generados:")
    for i in range(1, NUM_RUNS + 1):
        output_file = f"output_run{i}_scored.mol2"
        if os.path.exists(output_file):
            size = os.path.getsize(output_file)
            print(f"  ✓ {output_file} ({size} bytes)")
        else:
            print(f"  ✗ {output_file} - NO GENERADO")
    
    print("\n" + "=" * 70)
    print("Ejecuta 'python3 analyze_dock6_results.py' para calcular RMSD")
    print("=" * 70)


if __name__ == "__main__":
    main()
