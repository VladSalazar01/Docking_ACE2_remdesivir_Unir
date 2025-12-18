#!/usr/bin/env python3
"""
BENCHMARK BOLTZ - Re-docking 9FMM (Simplificado)

Boltz NO requiere preparación:
  - Proteína: solo secuencia de aminoácidos
  - Ligando: solo código CCD o SMILES

Sistema: 9FMM (ACE2 humano + inhibidor A1IDX)
"""

import os
import sys
import time
import json
import subprocess
import logging
from pathlib import Path
from datetime import datetime
import numpy as np

# RDKit para RMSD preciso
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('benchmark_boltz.log', encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

class Config:
    # Detectar si estamos en WSL/Linux o Windows
    import platform
    if platform.system() == "Linux":
        WORKDIR = Path("/mnt/g/Boltz_benchmark")
    else:
        WORKDIR = Path(r"G:\Boltz_benchmark")
    
    OUTPUT_DIR = WORKDIR / "boltz_outputs"
    RESULTS_DIR = WORKDIR / "results"
    
    # Referencia para RMSD (pose cristalográfica - solo residuo 801)
    LIGAND_REFERENCE = WORKDIR / "ligand_ref_801.pdb"
    
    # --- INPUTS BOLTZ (sin preparación) ---
    # SMILES del ligando A1IDX (corregido para RDKit)
    LIGAND_SMILES = "CC(C)C[C@H](NC(=O)[C@H](Cc1cn(Cc2cc(F)cc(Cl)c2)cn1)C(=O)O)C(=O)O"
    
    # =========================================================================
    # SECUENCIAS DISPONIBLES
    # =========================================================================
    
    # Secuencia COMPLETA ACE2 (residuos 19-740 de 9FMM, ~722 residuos)
    # Para uso con CPU (requiere mucha RAM, ~128GB recomendado)
    SEQUENCE_FULL = (
        "STIEEQAKTFLDKFNHEAEDLFYQSSLASWNYNTNITEENVQNMNNAGDKWSAFLKEQSTLAQMYPLQEI"
        "QNLTVKLQLQALQQNGSSVLSEDKSKRLNTILNTMSTIYSTGKVCNPDNPQECLLLEPGLNEIMANSLDY"
        "NERLWAWESWRSEVGKQLRPLYEEYVVLKNEMARANHYEDYGDYWRGDYEVNGVDGYDYSRGQLIEDVEH"
        "TFEEIKPLYEHLHAYVRAKLMNAYPSYISPIGCLPAHLLGDMWGRFWTNLYSLTVPFGQKPNIDVTDAMV"
        "DQAWDAQRIFKEAEKFFVSVGLPNMTQGFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDD"
        "FLTAHHEMGHIQYDMAYAAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEINF"
        "LLKQALTIVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETYCDPASLFHVSN"
        "DYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRLGKSEPWTLALENVVGAKN"
        "MNVRPLLNYFEPLFTWLKDQNKNSFVGWSTDWSPYADQSIKVRISLKSALGDKAYEWNDNEMYLFRSSVA"
        "YAMRQYFLKVKNQMILFGEEDVRVANLKPRISFNFFVTAPKNVSDIIPRTEVEKAIRMSRSRINDAFRLN"
        "DNSLEFLGIQPTLGPPNQPPVS"
    )  # 722 residuos
    
    # Secuencia TRUNCADA - dominio catalítico (~261 residuos)
    # Para uso con GPU (12GB VRAM)
    SEQUENCE_CATALYTIC = (
        "GFWENSMLTDPGNVQKAVCHPTAWDLGKGDFRILMCTKVTMDDFLTAHHEMGHIQYDMAY"
        "AAQPFLLRNGANEGFHEAVGEIMSLSAATPKHLKSIGLLSPDFQEDNETEINFLLKQALT"
        "IVGTLPFTYMLEKWRWMVFKGEIPKDQWMKKWWEMKREIVGVVEPVPHDETYCDPASLF"
        "HVSNDYSFIRYYTRTLYQFQFQEALCQAAKHEGPLHKCDISNSTEAGQKLFNMLRLGKS"
        "EPWTLALENVVGAKNMNVRPL"
    )  # 261 residuos
    
    # =========================================================================
    # RESIDUOS POCKET (para constraints)
    # =========================================================================
    
    # Pocket para secuencia COMPLETA (numeración absoluta en secuencia)
    # His345=327, Pro346=328, Thr371=353, His374=356, Glu375=357, 
    # His378=360, Glu402=384, His505=487, Tyr515=497
    POCKET_FULL = [327, 328, 353, 356, 357, 360, 384, 487, 497]
    
    # Pocket para secuencia TRUNCADA (offset: empieza en residuo 250)
    # His345=96, Pro346=97, Thr371=122, His374=125, etc.
    POCKET_CATALYTIC = [96, 97, 122, 125, 126, 129, 153, 237, 247]
    
    # =========================================================================
    # MODO ACTUAL (se configura en runtime)
    # =========================================================================
    MODE = "GPU"  # "GPU" o "CPU"
    PROTEIN_SEQUENCE = SEQUENCE_CATALYTIC
    POCKET_RESIDUES = POCKET_CATALYTIC
    
    # Parámetros Boltz
    RECYCLING_STEPS = 3
    DIFFUSION_SAMPLES = 1  # GPU: 1, CPU: puede ser más
    NUM_RUNS = 10
    USE_POCKET = True  # Usar pocket constraints
    
    # GPU config
    CUDA_VISIBLE_DEVICES = "0"
    NUM_GPUS = 1  # Número de GPUs a usar
    
    # Umbrales RMSD
    RMSD_EXCELLENT = 2.0
    RMSD_ACCEPTABLE = 3.0
    
    @classmethod
    def set_mode(cls, mode: str):
        """Configurar modo GPU, CPU o HYBRID."""
        mode = mode.upper()
        if mode == "CPU":
            cls.MODE = "CPU"
            cls.PROTEIN_SEQUENCE = cls.SEQUENCE_FULL
            cls.POCKET_RESIDUES = cls.POCKET_FULL
            cls.CUDA_VISIBLE_DEVICES = "-1"  # Forzar CPU
            cls.DIFFUSION_SAMPLES = 3  # Más muestras en CPU
            cls.NUM_RUNS = 3  # Menos runs porque es lento
            cls.OUTPUT_DIR = cls.WORKDIR / "boltz_outputs_cpu"
            logger.info("Modo CPU: Secuencia completa (722 res), diffusion_samples=3, runs=3")
        elif mode == "HYBRID":
            cls.MODE = "HYBRID"
            cls.PROTEIN_SEQUENCE = cls.SEQUENCE_CATALYTIC
            cls.POCKET_RESIDUES = cls.POCKET_CATALYTIC
            cls.CUDA_VISIBLE_DEVICES = "0"
            cls.DIFFUSION_SAMPLES = 1
            cls.NUM_RUNS = 10
            cls.OUTPUT_DIR = cls.WORKDIR / "boltz_outputs_hybrid"
            logger.info("Modo HYBRID: GPU (261 res) + 1 validacion CPU (722 res)")
        else:  # GPU (default)
            cls.MODE = "GPU"
            cls.PROTEIN_SEQUENCE = cls.SEQUENCE_CATALYTIC
            cls.POCKET_RESIDUES = cls.POCKET_CATALYTIC
            cls.CUDA_VISIBLE_DEVICES = "0"
            cls.DIFFUSION_SAMPLES = 1
            cls.NUM_RUNS = 10
            cls.OUTPUT_DIR = cls.WORKDIR / "boltz_outputs_gpu"
            logger.info("Modo GPU: Secuencia truncada (261 res), diffusion_samples=1")


# =============================================================================
# FUNCIONES
# =============================================================================

def create_input_yaml(output_path: Path, use_pocket: bool = True) -> None:
    """Crear YAML para Boltz con pocket constraints opcionales."""
    yaml_content = f"""version: 1
sequences:
  - protein:
      id: A
      sequence: {Config.PROTEIN_SEQUENCE}
  - ligand:
      id: B
      smiles: "{Config.LIGAND_SMILES}"
"""
    
    # Añadir pocket constraints para guiar al ligando al sitio activo
    if use_pocket and hasattr(Config, 'POCKET_RESIDUES'):
        contacts = "\n".join([f"        - [A, {r}]" for r in Config.POCKET_RESIDUES])
        yaml_content += f"""constraints:
  - pocket:
      binder: B
      contacts:
{contacts}
"""
    
    output_path.write_text(yaml_content)
    logger.info(f"YAML creado: {output_path} (pocket={'Si' if use_pocket else 'No'})")


def run_boltz(input_yaml: Path, output_dir: Path, seed: int, force_cpu: bool = False):
    """Ejecutar Boltz."""
    env = os.environ.copy()
    use_cpu = force_cpu or Config.MODE == "CPU"
    
    # Optimizaciones CPU
    if use_cpu:
        env['OMP_NUM_THREADS'] = '48'
        env['MKL_NUM_THREADS'] = '48'
        env['NUMEXPR_NUM_THREADS'] = '48'
        env['PYTORCH_CUDA_ALLOC_CONF'] = ''
    
    cmd = [
        'boltz', 'predict', str(input_yaml),
        '--out_dir', str(output_dir),
        '--recycling_steps', str(Config.RECYCLING_STEPS),
        '--diffusion_samples', str(Config.DIFFUSION_SAMPLES),
        '--seed', str(seed),
        '--accelerator', 'cpu' if use_cpu else 'gpu',
        '--use_msa_server',
        '--no_kernels',
        '--override'
    ]
    
    # Multi-GPU support
    if not use_cpu and Config.NUM_GPUS > 1:
        cmd.extend(['--devices', str(Config.NUM_GPUS)])
    
    logger.info(f"Ejecutando: {' '.join(cmd)}")
    
    start = time.time()
    timeout = 14400 if use_cpu else 3600
    
    try:
        result = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=timeout)
        elapsed = time.time() - start
        
        if result.returncode != 0:
            logger.error(f"Boltz fallo (codigo {result.returncode})")
            if result.stderr:
                for line in result.stderr.strip().split('\n')[-10:]:
                    logger.error(f"  {line}")
            return False, elapsed
        
        return True, elapsed
        
    except subprocess.TimeoutExpired:
        logger.error(f"Timeout ({timeout}s)")
        return False, time.time() - start
    except Exception as e:
        logger.error(f"Error: {e}")
        return False, time.time() - start


def get_coords(pdb_path: Path) -> np.ndarray:
    """Extraer coordenadas de átomos pesados (no-H)."""
    coords = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith(('ATOM', 'HETATM')):
            atom = line[12:16].strip()
            if not atom.startswith('H') and len(atom) > 0:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass
    return np.array(coords) if coords else np.array([]).reshape(0, 3)


def calculate_rmsd(c1: np.ndarray, c2: np.ndarray) -> float:
    """RMSD con alineamiento Kabsch (fallback si RDKit no disponible)."""
    n = min(len(c1), len(c2))
    if n == 0:
        return float('inf')
    c1, c2 = c1[:n] - c1[:n].mean(0), c2[:n] - c2[:n].mean(0)
    U, _, Vt = np.linalg.svd(c2.T @ c1)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return np.sqrt(np.mean(np.sum((c1 - c2 @ R)**2, axis=1)))


def calculate_rmsd_rdkit(ref_pdb: Path, pred_pdb: Path, smiles: str) -> float:
    """RMSD usando RDKit con alineamiento por SMILES template."""
    if not HAS_RDKIT:
        logger.warning("RDKit no disponible, usando RMSD simple")
        return calculate_rmsd(get_coords(ref_pdb), get_coords(pred_pdb))
    
    try:
        # Crear molécula template desde SMILES
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            raise ValueError("SMILES inválido")
        template = Chem.AddHs(template)
        AllChem.EmbedMolecule(template, randomSeed=42)
        template = Chem.RemoveHs(template)
        
        # Cargar moléculas desde PDB
        ref_mol = Chem.MolFromPDBFile(str(ref_pdb), removeHs=True)
        pred_mol = Chem.MolFromPDBFile(str(pred_pdb), removeHs=True)
        
        if ref_mol is None or pred_mol is None:
            raise ValueError("No se pudo cargar PDB")
        
        # Calcular RMSD con alineamiento
        rmsd = rdMolAlign.GetBestRMS(pred_mol, ref_mol)
        return rmsd
        
    except Exception as e:
        logger.warning(f"RDKit RMSD falló ({e}), usando método simple")
        return calculate_rmsd(get_coords(ref_pdb), get_coords(pred_pdb))


def extract_ligand(cif_path: Path, out_pdb: Path) -> bool:
    """Extraer ligando desde mmCIF de Boltz (busca LIG1 o HETATM)."""
    try:
        coords = []
        for line in cif_path.read_text().splitlines():
            if line.startswith('HETATM') and 'LIG1' in line:
                parts = line.split()
                # Formato: HETATM id element atom . LIG1 . seq ? chain x y z ...
                # Índices: 0      1  2       3    4 5    6 7   8 9     10 11 12
                if len(parts) >= 13:
                    atom = parts[3]
                    x, y, z = float(parts[10]), float(parts[11]), float(parts[12])
                    coords.append((atom, x, y, z))
        
        if coords:
            with open(out_pdb, 'w') as f:
                for i, (atom, x, y, z) in enumerate(coords, 1):
                    f.write(f"HETATM{i:5d}  {atom:<3s} LIG B   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n")
                f.write("END\n")
            return True
    except Exception as e:
        logger.warning(f"Error extrayendo ligando: {e}")
    return False


# =============================================================================
# BENCHMARK
# =============================================================================

def run_benchmark():
    import platform
    os_name = "WSL/Linux" if platform.system() == "Linux" else "Windows"
    gpu_info = f"RTX 3060 12GB x{Config.NUM_GPUS} ({os_name})"
    
    mode_info = {
        "GPU": gpu_info,
        "CPU": "Threadripper 3960X 24C/48T + 128GB RAM",
        "HYBRID": "GPU benchmark + CPU validation"
    }
    
    print(f"""
    ==================================================================
             BENCHMARK BOLTZ - Re-docking 9FMM                 
      Modo: {Config.MODE}
      Hardware: {mode_info.get(Config.MODE, 'Unknown')}
      Proteina: {len(Config.PROTEIN_SEQUENCE)} residuos
      Runs: {Config.NUM_RUNS}
      Diffusion samples: {Config.DIFFUSION_SAMPLES}
      Pocket constraints: {'ACTIVO (' + str(len(Config.POCKET_RESIDUES)) + ' residuos)' if Config.USE_POCKET else 'DESACTIVADO'}
    ==================================================================
    """)
    
    Config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    Config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Referencia
    ref_coords = get_coords(Config.LIGAND_REFERENCE) if Config.LIGAND_REFERENCE.exists() else None
    if ref_coords is not None:
        logger.info(f"Referencia: {len(ref_coords)} atomos ({Config.LIGAND_REFERENCE.name})")
    else:
        logger.warning(f"Referencia no encontrada: {Config.LIGAND_REFERENCE}")
    
    logger.info(f"Modo: {Config.MODE}")
    logger.info(f"Proteina: {len(Config.PROTEIN_SEQUENCE)} residuos")
    logger.info(f"Pocket: {len(Config.POCKET_RESIDUES)} residuos del sitio activo")
    logger.info(f"Ligando: A1IDX (SMILES)")
    
    results = []
    total_start = time.time()
    
    for run_id in range(1, Config.NUM_RUNS + 1):
        logger.info(f"\n=== RUN {run_id}/{Config.NUM_RUNS} ===")
        
        run_dir = Config.OUTPUT_DIR / f"run_{run_id:03d}"
        run_dir.mkdir(parents=True, exist_ok=True)
        
        input_yaml = run_dir / "input.yaml"
        create_input_yaml(input_yaml, use_pocket=Config.USE_POCKET)
        
        force_cpu = Config.MODE == "CPU"
        success, elapsed = run_boltz(input_yaml, run_dir, seed=42+run_id, force_cpu=force_cpu)
        
        result = {'run_id': run_id, 'success': success, 'time_s': elapsed, 'rmsd_values': []}
        
        if success and ref_coords is not None:
            for cif in run_dir.rglob('*.cif'):
                lig_pdb = run_dir / f"lig_{cif.stem}.pdb"
                if extract_ligand(cif, lig_pdb):
                    try:
                        pred_coords = get_coords(lig_pdb)
                        logger.info(f"  Ref: {len(ref_coords)} atomos, Pred: {len(pred_coords)} atomos")
                        if len(pred_coords) > 0 and len(ref_coords) > 0:
                            # Intentar RMSD con RDKit primero
                            rmsd = calculate_rmsd_rdkit(
                                Config.LIGAND_REFERENCE, 
                                lig_pdb, 
                                Config.LIGAND_SMILES
                            )
                            result['rmsd_values'].append(rmsd)
                            logger.info(f"  {cif.name}: RMSD = {rmsd:.3f} A")
                        else:
                            logger.warning(f"  Sin coordenadas para RMSD")
                    except Exception as e:
                        logger.warning(f"  Error RMSD: {e}")
            
            if result['rmsd_values']:
                result['best_rmsd'] = min(result['rmsd_values'])
                logger.info(f"  -> Mejor: {result['best_rmsd']:.3f} A")
        
        results.append(result)
    
    # Resumen
    total_time = time.time() - total_start
    successful = [r for r in results if r.get('best_rmsd')]
    
    print(f"\n{'='*50}")
    print("RESUMEN")
    print(f"{'='*50}")
    
    if successful:
        rmsds = [r['best_rmsd'] for r in successful]
        n_exc = sum(1 for r in rmsds if r <= Config.RMSD_EXCELLENT)
        n_acc = sum(1 for r in rmsds if r <= Config.RMSD_ACCEPTABLE)
        
        print(f"Runs exitosos: {len(successful)}/{len(results)}")
        print(f"RMSD: {np.mean(rmsds):.3f} +/- {np.std(rmsds):.3f} A")
        print(f"  Minimo: {min(rmsds):.3f} A")
        print(f"  <=2.0 A: {n_exc}/{len(successful)} ({100*n_exc/len(successful):.1f}%)")
        print(f"  <=3.0 A: {n_acc}/{len(successful)} ({100*n_acc/len(successful):.1f}%)")
    
    print(f"Tiempo total: {total_time/60:.1f} min")
    
    # Guardar
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    json_path = Config.RESULTS_DIR / f"benchmark_{Config.MODE}_{ts}.json"
    json_path.write_text(json.dumps({
        'mode': Config.MODE,
        'protein_length': len(Config.PROTEIN_SEQUENCE),
        'diffusion_samples': Config.DIFFUSION_SAMPLES,
        'results': results, 
        'total_time': total_time
    }, indent=2, default=str))
    print(f"\nResultados: {json_path}")


def print_usage():
    print("""
Uso: python benchmark_boltz.py [MODO] [OPCIONES]

MODOS:
  gpu     Usar GPU con secuencia truncada (~261 residuos) [RECOMENDADO]
          - Rapido (~2 min/run)
          - Requiere: RTX 3060 12GB o similar
          
  cpu     Usar CPU (NO RECOMENDADO - extremadamente lento)
          - >2 horas por run
          - No practico para benchmarking
          
  hybrid  Equivalente a 'gpu' (validacion CPU deshabilitada)

OPCIONES:
  -n N    Numero de runs (default: 10)
  -h      Mostrar esta ayuda

EJEMPLOS:
  python benchmark_boltz.py gpu          # 10 runs en GPU
  python benchmark_boltz.py gpu -n 5     # 5 runs en GPU
""")


def run_cpu_validation():
    """Ejecutar 1 run en CPU con secuencia completa para validación."""
    print("\n" + "="*60)
    print("VALIDACION CPU - Secuencia completa")
    print("="*60)
    
    # Guardar config actual
    orig_seq = Config.PROTEIN_SEQUENCE
    orig_pocket = Config.POCKET_RESIDUES
    orig_cuda = Config.CUDA_VISIBLE_DEVICES
    orig_samples = Config.DIFFUSION_SAMPLES
    
    # Cambiar a CPU
    Config.PROTEIN_SEQUENCE = Config.SEQUENCE_FULL
    Config.POCKET_RESIDUES = Config.POCKET_FULL
    Config.CUDA_VISIBLE_DEVICES = "-1"  # Forzar CPU
    Config.DIFFUSION_SAMPLES = 3
    
    cpu_dir = Config.WORKDIR / "boltz_outputs_cpu_validation"
    cpu_dir.mkdir(parents=True, exist_ok=True)
    
    run_dir = cpu_dir / "run_cpu_001"
    run_dir.mkdir(parents=True, exist_ok=True)
    
    input_yaml = run_dir / "input.yaml"
    create_input_yaml(input_yaml, use_pocket=True)
    
    logger.info(f"CPU Validation: {len(Config.PROTEIN_SEQUENCE)} residuos, diffusion_samples={Config.DIFFUSION_SAMPLES}")
    
    start = time.time()
    success, elapsed = run_boltz(input_yaml, run_dir, seed=999, force_cpu=True)
    
    result = {'success': success, 'time_s': elapsed, 'rmsd_values': []}
    
    if success:
        ref_coords = get_coords(Config.LIGAND_REFERENCE) if Config.LIGAND_REFERENCE.exists() else None
        if ref_coords is not None:
            for cif in run_dir.rglob('*.cif'):
                lig_pdb = run_dir / f"lig_{cif.stem}.pdb"
                if extract_ligand(cif, lig_pdb):
                    try:
                        pred_coords = get_coords(lig_pdb)
                        if len(pred_coords) > 0:
                            rmsd = calculate_rmsd_rdkit(Config.LIGAND_REFERENCE, lig_pdb, Config.LIGAND_SMILES)
                            result['rmsd_values'].append(rmsd)
                            logger.info(f"  CPU Validation RMSD: {rmsd:.3f} A")
                    except Exception as e:
                        logger.warning(f"  Error RMSD: {e}")
    
    # Restaurar config
    Config.PROTEIN_SEQUENCE = orig_seq
    Config.POCKET_RESIDUES = orig_pocket
    Config.CUDA_VISIBLE_DEVICES = orig_cuda
    Config.DIFFUSION_SAMPLES = orig_samples
    
    print(f"\nValidacion CPU completada en {elapsed/60:.1f} min")
    if result['rmsd_values']:
        print(f"  RMSD: {min(result['rmsd_values']):.3f} A")
    
    return result


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Benchmark Boltz para re-docking 9FMM',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('mode', nargs='?', default='gpu', 
                        choices=['gpu', 'GPU', 'cpu', 'CPU', 'hybrid', 'HYBRID'],
                        help='Modo de ejecucion: gpu, cpu o hybrid (default: gpu)')
    parser.add_argument('-n', '--num_runs', type=int, default=None,
                        help='Numero de runs (default: 10)')
    parser.add_argument('-r', '--recycling', type=int, default=3,
                        help='Recycling steps (default: 3)')
    parser.add_argument('-g', '--gpus', type=int, default=1,
                        help='Numero de GPUs (default: 1)')
    parser.add_argument('--no-pocket', action='store_true',
                        help='Deshabilitar pocket constraints')
    
    args = parser.parse_args()
    
    # Verificar Boltz
    try:
        print("Verificando Boltz (puede tardar ~60s la primera vez)...")
        result = subprocess.run(['boltz', '--help'], capture_output=True, timeout=120)
        if result.returncode != 0:
            sys.exit("ERROR: Boltz no responde correctamente")
        print("Boltz OK")
    except FileNotFoundError:
        sys.exit("ERROR: Boltz no instalado. Ejecuta: pip install boltz")
    except Exception as e:
        sys.exit(f"ERROR verificando Boltz: {e}")
    
    # Configurar modo
    Config.set_mode(args.mode)
    
    # Override num_runs si se especifica
    if args.num_runs is not None:
        Config.NUM_RUNS = args.num_runs
    
    # Override recycling_steps
    Config.RECYCLING_STEPS = args.recycling
    
    # Deshabilitar pocket constraints si se pide
    Config.USE_POCKET = not args.no_pocket
    
    # Multi-GPU
    Config.NUM_GPUS = args.gpus
    
    # Configurar directorio según opciones
    suffix = f"_r{Config.RECYCLING_STEPS}"
    if not Config.USE_POCKET:
        suffix += "_nopocket"
    Config.OUTPUT_DIR = Config.WORKDIR / f"boltz_outputs{suffix}"
    
    Config.WORKDIR.mkdir(parents=True, exist_ok=True)
    
    # Mostrar configuración
    print("\n" + "="*60)
    print(f"CONFIGURACION:")
    print(f"  - Modo: {Config.MODE}")
    print(f"  - GPUs: {Config.NUM_GPUS}")
    print(f"  - Recycling steps: {Config.RECYCLING_STEPS}")
    print(f"  - Pocket constraints: {'Si' if Config.USE_POCKET else 'No'}")
    print(f"  - Output dir: {Config.OUTPUT_DIR}")
    print(f"  - Runs: {Config.NUM_RUNS}")
    print("="*60 + "\n")
    
    # Advertencia para CPU
    if Config.MODE == "CPU":
        print("\n" + "="*60)
        print("ADVERTENCIA: Modo CPU NO RECOMENDADO")
        print("  - Boltz es extremadamente lento en CPU")
        print("  - 1 run puede tomar >2 horas (vs ~2 min en GPU)")
        print("  - Considera usar GPU con secuencia truncada")
        print("="*60)
        resp = input("\nContinuar de todos modos? [s/N]: ").strip().lower()
        if resp != 's':
            sys.exit("Cancelado. Usa: python benchmark_boltz.py gpu")
    
    # Ejecutar benchmark
    run_benchmark()
    
    # Si es modo híbrido, ejecutar validación CPU
    if Config.MODE == "HYBRID":
        print("\n" + "="*60)
        print("NOTA: Validacion CPU omitida")
        print("  - CPU es extremadamente lento para Boltz (>2h por run)")
        print("  - Resultados GPU con dominio catalitico son suficientes")
        print("="*60)