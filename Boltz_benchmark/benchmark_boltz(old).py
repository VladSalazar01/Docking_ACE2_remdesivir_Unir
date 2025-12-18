#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
BENCHMARK BOLTZ - Re-docking Pipeline para 9FMM
=============================================================================

Script de benchmarking para Boltz-2 (predicción de estructuras proteína-ligando).
Adaptado del pipeline de benchmark_gnina_dual_site.py para el TFM.

Sistema: 9FMM (ACE2 humano + inhibidor fluorinado A1IDX)
Hardware: AMD Ryzen Threadripper 3690x 24c/48t, 128GB RAM, 2x RTX3060

Autor: Adaptación para TFM
Fecha: 2025

Requisitos:
    - boltz (pip install boltz[cuda])
    - rdkit
    - numpy
    - biopython
    
Uso:
    python benchmark_boltz.py
"""

import os
import sys
import time
import json
import yaml
import shutil
import subprocess
import logging
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import numpy as np

# Configuración de logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('benchmark_boltz.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# =============================================================================
# CONFIGURACIÓN
# =============================================================================

class Config:
    """Configuración del benchmark para Boltz."""
    
    # Directorio de trabajo (Windows path)
    WORKDIR = Path(r"G:\Boltz_benchmark")
    
    # Archivos de entrada
    PROTEIN_PDB = "9FMM_protein.pdb"
    LIGAND_PDB = "ligand_A1IDX.pdb"
    CIF_FILE = "9FMM.cif"
    
    # Código CCD del ligando (desde PDB)
    LIGAND_CCD = "A1IDX"
    
    # SMILES del ligando A1IDX (F-MLN-4760)
    # Estructura: (2S)-2-[[(2S)-3-[3-[(3-chloranyl-5-fluoranyl-phenyl)methyl]imidazol-4-yl]-1-oxidanyl-1-oxidanylidene-propan-2-yl]amino]-4-methyl-pentanoic acid
    LIGAND_SMILES = "CC(C)C[C@H](NC(=O)[C@H](Cc1c[nH]cn1Cc2cc(F)cc(Cl)c2)C(=O)O)C(=O)O"
    
    # Parámetros de Boltz
    RECYCLING_STEPS = 3        # Pasos de refinamiento (default: 3)
    DIFFUSION_SAMPLES = 5      # Número de muestras/predicciones (para ensemble)
    SAMPLING_STEPS = 50        # Pasos de difusión (default: 50)
    
    # Número de runs de benchmark
    NUM_RUNS = 10
    
    # Directorio de salida
    OUTPUT_DIR = WORKDIR / "boltz_outputs"
    RESULTS_DIR = WORKDIR / "benchmark_results"
    
    # Threshold para éxito de RMSD
    RMSD_SUCCESS_THRESHOLD = 2.0  # Å
    RMSD_ACCEPTABLE_THRESHOLD = 3.0  # Å
    
    # GPU configuration
    CUDA_VISIBLE_DEVICES = "0,1"  # Usar ambas RTX3060


# =============================================================================
# FUNCIONES AUXILIARES
# =============================================================================

def extract_sequence_from_pdb(pdb_path: Path) -> str:
    """
    Extraer la secuencia de aminoácidos desde un archivo PDB.
    
    Args:
        pdb_path: Ruta al archivo PDB
        
    Returns:
        Secuencia de aminoácidos en formato de una letra
    """
    try:
        from Bio.PDB import PDBParser
        from Bio.SeqUtils import seq1
    except ImportError:
        logger.error("BioPython no instalado. Ejecutar: pip install biopython")
        raise
    
    # Código de tres letras a una letra
    aa_map = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    
    sequences = []
    for model in structure:
        for chain in model:
            seq = ""
            for residue in chain:
                res_name = residue.get_resname().strip()
                if res_name in aa_map:
                    seq += aa_map[res_name]
            if seq:
                sequences.append(seq)
    
    # Retornar la secuencia más larga (cadena principal)
    if sequences:
        return max(sequences, key=len)
    else:
        raise ValueError(f"No se encontró secuencia en {pdb_path}")


def extract_ligand_smiles_from_pdb(pdb_path: Path, ligand_name: str = "A1IDX") -> str:
    """
    Extraer SMILES del ligando desde un archivo PDB usando RDKit.
    
    Args:
        pdb_path: Ruta al archivo PDB del ligando
        ligand_name: Nombre del ligando
        
    Returns:
        SMILES string
    """
    try:
        from rdkit import Chem
    except ImportError:
        logger.warning("RDKit no disponible. Usando SMILES predefinido.")
        return Config.LIGAND_SMILES
    
    try:
        mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False)
        if mol:
            smiles = Chem.MolToSmiles(mol)
            logger.info(f"SMILES extraído del ligando: {smiles}")
            return smiles
    except Exception as e:
        logger.warning(f"Error extrayendo SMILES: {e}. Usando predefinido.")
    
    return Config.LIGAND_SMILES


def get_ligand_coords_from_pdb(pdb_path: Path) -> np.ndarray:
    """
    Obtener coordenadas del ligando desde archivo PDB.
    
    Args:
        pdb_path: Ruta al archivo PDB del ligando
        
    Returns:
        Array numpy con coordenadas (N, 3)
    """
    coords = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extraer solo átomos pesados (no hidrógenos)
                atom_name = line[12:16].strip()
                if not atom_name.startswith('H'):
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
    
    return np.array(coords)


def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Calcular RMSD entre dos conjuntos de coordenadas usando alineamiento de Kabsch.
    
    Args:
        coords1: Coordenadas de referencia (N, 3)
        coords2: Coordenadas predichas (N, 3)
        
    Returns:
        RMSD en Angstroms
    """
    # Asegurar mismo número de átomos (usar el mínimo)
    n = min(len(coords1), len(coords2))
    c1 = coords1[:n]
    c2 = coords2[:n]
    
    # Centrar coordenadas
    c1_centered = c1 - c1.mean(axis=0)
    c2_centered = c2 - c2.mean(axis=0)
    
    # Matriz de covarianza
    H = c2_centered.T @ c1_centered
    
    # SVD para rotación óptima
    U, S, Vt = np.linalg.svd(H)
    
    # Matriz de rotación (con corrección de reflexión)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    
    # Aplicar rotación
    c2_aligned = c2_centered @ R
    
    # Calcular RMSD
    diff = c1_centered - c2_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    
    return rmsd


def calculate_rmsd_mcs(ref_pdb: Path, pred_pdb: Path) -> float:
    """
    Calcular RMSD usando Maximum Common Substructure (MCS).
    Método topology-aware para moléculas con simetría.
    
    Args:
        ref_pdb: Ruta al PDB de referencia
        pred_pdb: Ruta al PDB predicho
        
    Returns:
        RMSD en Angstroms
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign
        
        ref_mol = Chem.MolFromPDBFile(str(ref_pdb), removeHs=True)
        pred_mol = Chem.MolFromPDBFile(str(pred_pdb), removeHs=True)
        
        if ref_mol and pred_mol:
            # Alineamiento con MCS
            rmsd = rdMolAlign.GetBestRMS(ref_mol, pred_mol)
            return rmsd
    except Exception as e:
        logger.warning(f"Error en cálculo MCS-RMSD: {e}. Usando método simple.")
    
    # Fallback a método simple
    ref_coords = get_ligand_coords_from_pdb(ref_pdb)
    pred_coords = get_ligand_coords_from_pdb(pred_pdb)
    return calculate_rmsd(ref_coords, pred_coords)


# =============================================================================
# GENERACIÓN DE INPUTS PARA BOLTZ
# =============================================================================

def create_boltz_yaml(
    output_path: Path,
    protein_sequence: str,
    ligand_ccd: str = None,
    ligand_smiles: str = None,
    protein_id: str = "A",
    ligand_id: str = "B"
) -> None:
    """
    Crear archivo YAML de entrada para Boltz.
    
    Args:
        output_path: Ruta donde guardar el archivo YAML
        protein_sequence: Secuencia de la proteína
        ligand_ccd: Código CCD del ligando (opcional)
        ligand_smiles: SMILES del ligando (opcional, alternativo a CCD)
        protein_id: ID de la cadena de proteína
        ligand_id: ID de la cadena del ligando
    """
    yaml_content = {
        'version': 1,
        'sequences': [
            {
                'protein': {
                    'id': protein_id,
                    'sequence': protein_sequence
                }
            }
        ]
    }
    
    # Agregar ligando (preferir CCD sobre SMILES)
    if ligand_ccd:
        yaml_content['sequences'].append({
            'ligand': {
                'id': ligand_id,
                'ccd': ligand_ccd
            }
        })
        logger.info(f"Usando ligando CCD: {ligand_ccd}")
    elif ligand_smiles:
        yaml_content['sequences'].append({
            'ligand': {
                'id': ligand_id,
                'smiles': ligand_smiles
            }
        })
        logger.info(f"Usando ligando SMILES: {ligand_smiles[:50]}...")
    else:
        raise ValueError("Debe proporcionar ligand_ccd o ligand_smiles")
    
    # Escribir YAML
    with open(output_path, 'w') as f:
        yaml.dump(yaml_content, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Archivo YAML creado: {output_path}")


# =============================================================================
# EJECUCIÓN DE BOLTZ
# =============================================================================

def run_boltz_prediction(
    input_yaml: Path,
    output_dir: Path,
    recycling_steps: int = 3,
    diffusion_samples: int = 5,
    sampling_steps: int = 50,
    use_msa_server: bool = True,
    seed: int = None
) -> Tuple[bool, float, dict]:
    """
    Ejecutar predicción con Boltz.
    
    Args:
        input_yaml: Ruta al archivo YAML de entrada
        output_dir: Directorio de salida
        recycling_steps: Número de pasos de refinamiento
        diffusion_samples: Número de muestras de difusión
        sampling_steps: Número de pasos de sampling
        use_msa_server: Usar servidor MSA (ColabFold)
        seed: Semilla aleatoria (opcional)
        
    Returns:
        Tuple (éxito, tiempo_ejecución, métricas)
    """
    # Configurar entorno
    env = os.environ.copy()
    env['CUDA_VISIBLE_DEVICES'] = Config.CUDA_VISIBLE_DEVICES
    env['TORCH_FORCE_FLOAT32_MATMUL_PRECISION'] = 'medium'
    
    # Construir comando
    cmd = [
        'boltz', 'predict',
        str(input_yaml),
        '--out_dir', str(output_dir),
        '--recycling_steps', str(recycling_steps),
        '--diffusion_samples', str(diffusion_samples),
        '--sampling_steps', str(sampling_steps),
    ]
    
    if use_msa_server:
        cmd.append('--use_msa_server')
    
    if seed is not None:
        cmd.extend(['--seed', str(seed)])
    
    logger.info(f"Ejecutando comando: {' '.join(cmd)}")
    
    start_time = time.time()
    success = False
    metrics = {}
    
    try:
        result = subprocess.run(
            cmd,
            env=env,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hora timeout
        )
        
        elapsed_time = time.time() - start_time
        
        if result.returncode == 0:
            success = True
            logger.info(f"Predicción completada en {elapsed_time:.2f}s")
            
            # Parsear métricas de confianza si están disponibles
            metrics = parse_boltz_outputs(output_dir)
        else:
            logger.error(f"Error en Boltz: {result.stderr}")
            
    except subprocess.TimeoutExpired:
        elapsed_time = time.time() - start_time
        logger.error(f"Timeout después de {elapsed_time:.2f}s")
    except Exception as e:
        elapsed_time = time.time() - start_time
        logger.error(f"Error ejecutando Boltz: {e}")
    
    return success, elapsed_time, metrics


def parse_boltz_outputs(output_dir: Path) -> dict:
    """
    Parsear archivos de salida de Boltz para extraer métricas.
    
    Args:
        output_dir: Directorio con salidas de Boltz
        
    Returns:
        Diccionario con métricas (pTM, ipTM, pLDDT, etc.)
    """
    metrics = {}
    
    # Buscar archivos de confianza JSON
    for json_file in output_dir.rglob("*confidence*.json"):
        try:
            with open(json_file, 'r') as f:
                data = json.load(f)
                metrics.update(data)
        except Exception as e:
            logger.warning(f"Error leyendo {json_file}: {e}")
    
    # Buscar archivos de ranking
    for ranking_file in output_dir.rglob("*ranking*.json"):
        try:
            with open(ranking_file, 'r') as f:
                data = json.load(f)
                metrics['ranking'] = data
        except Exception as e:
            logger.warning(f"Error leyendo {ranking_file}: {e}")
    
    return metrics


def extract_ligand_from_boltz_output(
    output_dir: Path,
    ligand_chain: str = "B",
    sample_idx: int = 0
) -> Optional[Path]:
    """
    Extraer coordenadas del ligando desde salida de Boltz.
    
    Args:
        output_dir: Directorio de salida de Boltz
        ligand_chain: ID de la cadena del ligando
        sample_idx: Índice de la muestra a extraer
        
    Returns:
        Ruta al archivo PDB del ligando extraído
    """
    # Buscar archivos de predicción (mmCIF o PDB)
    pred_files = list(output_dir.rglob("*.cif")) + list(output_dir.rglob("*.pdb"))
    
    if not pred_files:
        logger.error(f"No se encontraron archivos de predicción en {output_dir}")
        return None
    
    # Ordenar por nombre para obtener consistencia
    pred_files = sorted(pred_files)
    
    if sample_idx >= len(pred_files):
        logger.warning(f"sample_idx {sample_idx} fuera de rango. Usando 0.")
        sample_idx = 0
    
    pred_file = pred_files[sample_idx]
    logger.info(f"Procesando predicción: {pred_file}")
    
    # Extraer ligando
    output_ligand = output_dir / f"ligand_predicted_{sample_idx}.pdb"
    
    try:
        if pred_file.suffix == '.cif':
            # Extraer desde mmCIF
            extract_ligand_from_cif(pred_file, output_ligand, ligand_chain)
        else:
            # Extraer desde PDB
            extract_ligand_from_pdb_chain(pred_file, output_ligand, ligand_chain)
        
        return output_ligand
    except Exception as e:
        logger.error(f"Error extrayendo ligando: {e}")
        return None


def extract_ligand_from_cif(cif_path: Path, output_pdb: Path, chain_id: str) -> None:
    """
    Extraer ligando desde archivo mmCIF.
    
    Args:
        cif_path: Ruta al archivo mmCIF
        output_pdb: Ruta al archivo PDB de salida
        chain_id: ID de la cadena del ligando
    """
    try:
        from Bio.PDB import MMCIFParser, PDBIO, Select
        
        class LigandSelect(Select):
            def accept_chain(self, chain):
                return chain.id == chain_id
            def accept_residue(self, residue):
                return True  # Aceptar todos los residuos de la cadena
        
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("pred", str(cif_path))
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_pdb), LigandSelect())
        
        logger.info(f"Ligando extraído a: {output_pdb}")
        
    except ImportError:
        # Fallback: extracción manual
        logger.warning("BioPython no disponible. Usando extracción manual.")
        extract_ligand_manual_cif(cif_path, output_pdb, chain_id)


def extract_ligand_manual_cif(cif_path: Path, output_pdb: Path, chain_id: str) -> None:
    """
    Extracción manual de ligando desde mmCIF (fallback).
    """
    atoms = []
    
    with open(cif_path, 'r') as f:
        in_atom_site = False
        columns = {}
        
        for line in f:
            if line.startswith('_atom_site.'):
                in_atom_site = True
                col_name = line.split('.')[1].strip()
                columns[col_name] = len(columns)
            elif in_atom_site and line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) > max(columns.values()):
                    chain = parts[columns.get('auth_asym_id', columns.get('label_asym_id', 0))]
                    if chain == chain_id:
                        atoms.append(line)
    
    # Escribir PDB
    with open(output_pdb, 'w') as f:
        for i, atom in enumerate(atoms):
            parts = atom.split()
            # Formato PDB simplificado
            f.write(f"HETATM{i+1:5d}  {parts[3]:4s} LIG {chain_id}{1:4d}    "
                   f"{float(parts[10]):8.3f}{float(parts[11]):8.3f}{float(parts[12]):8.3f}\n")
        f.write("END\n")


def extract_ligand_from_pdb_chain(pdb_path: Path, output_pdb: Path, chain_id: str) -> None:
    """
    Extraer ligando desde archivo PDB por cadena.
    """
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    with open(output_pdb, 'w') as f:
        for line in lines:
            if (line.startswith('ATOM') or line.startswith('HETATM')):
                if len(line) > 21 and line[21] == chain_id:
                    f.write(line)
        f.write("END\n")


# =============================================================================
# BENCHMARK PRINCIPAL
# =============================================================================

class BoltzBenchmark:
    """Clase principal para ejecutar benchmark de Boltz."""
    
    def __init__(self):
        """Inicializar benchmark."""
        self.config = Config()
        self.results = []
        
        # Crear directorios
        self.config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        self.config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
        
        # Cargar datos de referencia
        self._load_reference_data()
    
    def _load_reference_data(self):
        """Cargar datos de referencia (secuencia y ligando experimental)."""
        protein_path = self.config.WORKDIR / self.config.PROTEIN_PDB
        ligand_path = self.config.WORKDIR / self.config.LIGAND_PDB
        
        # Verificar archivos
        if not protein_path.exists():
            raise FileNotFoundError(f"No se encontró: {protein_path}")
        if not ligand_path.exists():
            raise FileNotFoundError(f"No se encontró: {ligand_path}")
        
        # Extraer secuencia de proteína
        logger.info("Extrayendo secuencia de proteína...")
        self.protein_sequence = extract_sequence_from_pdb(protein_path)
        logger.info(f"Secuencia extraída: {len(self.protein_sequence)} residuos")
        
        # Guardar coordenadas de referencia del ligando
        self.ref_ligand_path = ligand_path
        self.ref_ligand_coords = get_ligand_coords_from_pdb(ligand_path)
        logger.info(f"Ligando de referencia: {len(self.ref_ligand_coords)} átomos pesados")
    
    def run_single_benchmark(
        self,
        run_id: int,
        use_ccd: bool = True
    ) -> dict:
        """
        Ejecutar una sola iteración de benchmark.
        
        Args:
            run_id: ID del run
            use_ccd: Usar código CCD (True) o SMILES (False)
            
        Returns:
            Diccionario con resultados del run
        """
        logger.info(f"\n{'='*60}")
        logger.info(f"BENCHMARK RUN {run_id}")
        logger.info(f"{'='*60}")
        
        run_dir = self.config.OUTPUT_DIR / f"run_{run_id:03d}"
        run_dir.mkdir(parents=True, exist_ok=True)
        
        # Crear archivo YAML de entrada
        input_yaml = run_dir / "input.yaml"
        create_boltz_yaml(
            output_path=input_yaml,
            protein_sequence=self.protein_sequence,
            ligand_ccd=self.config.LIGAND_CCD if use_ccd else None,
            ligand_smiles=self.config.LIGAND_SMILES if not use_ccd else None,
            protein_id="A",
            ligand_id="B"
        )
        
        # Ejecutar predicción
        seed = 42 + run_id  # Semilla diferente por run
        success, elapsed_time, metrics = run_boltz_prediction(
            input_yaml=input_yaml,
            output_dir=run_dir,
            recycling_steps=self.config.RECYCLING_STEPS,
            diffusion_samples=self.config.DIFFUSION_SAMPLES,
            sampling_steps=self.config.SAMPLING_STEPS,
            use_msa_server=True,
            seed=seed
        )
        
        result = {
            'run_id': run_id,
            'success': success,
            'elapsed_time': elapsed_time,
            'seed': seed,
            'metrics': metrics,
            'rmsd_values': [],
            'best_rmsd': None
        }
        
        if success:
            # Calcular RMSD para cada muestra de difusión
            for sample_idx in range(self.config.DIFFUSION_SAMPLES):
                pred_ligand = extract_ligand_from_boltz_output(
                    run_dir, 
                    ligand_chain="B",
                    sample_idx=sample_idx
                )
                
                if pred_ligand and pred_ligand.exists():
                    try:
                        rmsd = calculate_rmsd_mcs(
                            self.ref_ligand_path,
                            pred_ligand
                        )
                        result['rmsd_values'].append(rmsd)
                        logger.info(f"  Sample {sample_idx}: RMSD = {rmsd:.3f} Å")
                    except Exception as e:
                        logger.error(f"  Error calculando RMSD sample {sample_idx}: {e}")
            
            if result['rmsd_values']:
                result['best_rmsd'] = min(result['rmsd_values'])
                result['mean_rmsd'] = np.mean(result['rmsd_values'])
                result['std_rmsd'] = np.std(result['rmsd_values'])
                logger.info(f"  Mejor RMSD: {result['best_rmsd']:.3f} Å")
                logger.info(f"  RMSD medio: {result['mean_rmsd']:.3f} ± {result['std_rmsd']:.3f} Å")
        
        return result
    
    def run_benchmark(self) -> None:
        """Ejecutar benchmark completo."""
        logger.info("\n" + "="*70)
        logger.info("INICIANDO BENCHMARK BOLTZ PARA 9FMM")
        logger.info("="*70)
        logger.info(f"Proteína: ACE2 humano ({len(self.protein_sequence)} residuos)")
        logger.info(f"Ligando: {self.config.LIGAND_CCD}")
        logger.info(f"Número de runs: {self.config.NUM_RUNS}")
        logger.info(f"Muestras por run: {self.config.DIFFUSION_SAMPLES}")
        logger.info(f"Recycling steps: {self.config.RECYCLING_STEPS}")
        logger.info(f"Sampling steps: {self.config.SAMPLING_STEPS}")
        logger.info(f"GPUs: {self.config.CUDA_VISIBLE_DEVICES}")
        logger.info("="*70 + "\n")
        
        start_time = time.time()
        
        for run_id in range(1, self.config.NUM_RUNS + 1):
            result = self.run_single_benchmark(run_id)
            self.results.append(result)
        
        total_time = time.time() - start_time
        
        # Generar reporte
        self._generate_report(total_time)
    
    def _generate_report(self, total_time: float) -> None:
        """Generar reporte de resultados."""
        logger.info("\n" + "="*70)
        logger.info("RESUMEN DE RESULTADOS")
        logger.info("="*70)
        
        successful_runs = [r for r in self.results if r['success'] and r['best_rmsd'] is not None]
        
        if not successful_runs:
            logger.error("No se completaron runs exitosos.")
            return
        
        best_rmsds = [r['best_rmsd'] for r in successful_runs]
        all_rmsds = []
        for r in successful_runs:
            all_rmsds.extend(r['rmsd_values'])
        
        times = [r['elapsed_time'] for r in self.results]
        
        # Estadísticas
        stats = {
            'n_runs': len(self.results),
            'n_successful': len(successful_runs),
            'success_rate': len(successful_runs) / len(self.results) * 100,
            'best_rmsd_overall': min(best_rmsds),
            'worst_rmsd_overall': max(best_rmsds),
            'mean_best_rmsd': np.mean(best_rmsds),
            'std_best_rmsd': np.std(best_rmsds),
            'mean_all_rmsd': np.mean(all_rmsds),
            'std_all_rmsd': np.std(all_rmsds),
            'n_success_2A': sum(1 for r in best_rmsds if r <= self.config.RMSD_SUCCESS_THRESHOLD),
            'n_success_3A': sum(1 for r in best_rmsds if r <= self.config.RMSD_ACCEPTABLE_THRESHOLD),
            'mean_time': np.mean(times),
            'total_time': total_time
        }
        
        # Porcentajes de éxito
        stats['pct_success_2A'] = stats['n_success_2A'] / len(successful_runs) * 100
        stats['pct_success_3A'] = stats['n_success_3A'] / len(successful_runs) * 100
        
        # Imprimir resultados
        logger.info(f"Runs totales: {stats['n_runs']}")
        logger.info(f"Runs exitosos: {stats['n_successful']} ({stats['success_rate']:.1f}%)")
        logger.info("")
        logger.info(f"RMSD (mejor por run):")
        logger.info(f"  Mejor:   {stats['best_rmsd_overall']:.3f} Å")
        logger.info(f"  Peor:    {stats['worst_rmsd_overall']:.3f} Å")
        logger.info(f"  Media:   {stats['mean_best_rmsd']:.3f} ± {stats['std_best_rmsd']:.3f} Å")
        logger.info("")
        logger.info(f"RMSD (todas las muestras):")
        logger.info(f"  Media:   {stats['mean_all_rmsd']:.3f} ± {stats['std_all_rmsd']:.3f} Å")
        logger.info("")
        logger.info(f"Tasa de éxito:")
        logger.info(f"  RMSD ≤ 2.0 Å: {stats['n_success_2A']}/{len(successful_runs)} ({stats['pct_success_2A']:.1f}%)")
        logger.info(f"  RMSD ≤ 3.0 Å: {stats['n_success_3A']}/{len(successful_runs)} ({stats['pct_success_3A']:.1f}%)")
        logger.info("")
        logger.info(f"Tiempo:")
        logger.info(f"  Por run: {stats['mean_time']:.1f} s")
        logger.info(f"  Total:   {stats['total_time']:.1f} s ({stats['total_time']/60:.1f} min)")
        logger.info("="*70)
        
        # Guardar resultados en JSON
        results_path = self.config.RESULTS_DIR / f"benchmark_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        with open(results_path, 'w') as f:
            json.dump({
                'config': {
                    'protein': '9FMM (ACE2)',
                    'ligand': self.config.LIGAND_CCD,
                    'num_runs': self.config.NUM_RUNS,
                    'diffusion_samples': self.config.DIFFUSION_SAMPLES,
                    'recycling_steps': self.config.RECYCLING_STEPS,
                    'sampling_steps': self.config.SAMPLING_STEPS
                },
                'statistics': stats,
                'runs': self.results
            }, f, indent=2, default=str)
        
        logger.info(f"\nResultados guardados en: {results_path}")
        
        # Generar tabla CSV
        csv_path = self.config.RESULTS_DIR / f"benchmark_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        with open(csv_path, 'w') as f:
            f.write("Run,Success,Time(s),Best_RMSD(A),Mean_RMSD(A),Std_RMSD(A),N_Samples\n")
            for r in self.results:
                best = r['best_rmsd'] if r['best_rmsd'] else "N/A"
                mean = r.get('mean_rmsd', "N/A")
                std = r.get('std_rmsd', "N/A")
                n_samples = len(r['rmsd_values'])
                f.write(f"{r['run_id']},{r['success']},{r['elapsed_time']:.1f},"
                       f"{best},{mean},{std},{n_samples}\n")
        
        logger.info(f"Tabla CSV guardada en: {csv_path}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Función principal."""
    print("""
    ╔═══════════════════════════════════════════════════════════════╗
    ║           BENCHMARK BOLTZ - Re-docking 9FMM                   ║
    ║                                                               ║
    ║  Sistema: ACE2 humano + inhibidor fluorinado (A1IDX)          ║
    ║  Hardware: 2x RTX3060 + Threadripper 3690x                    ║
    ╚═══════════════════════════════════════════════════════════════╝
    """)
    
    try:
        # Verificar instalación de Boltz
        result = subprocess.run(['boltz', '--help'], capture_output=True, text=True)
        if result.returncode != 0:
            logger.error("Boltz no está instalado. Ejecutar: pip install boltz[cuda]")
            sys.exit(1)
        logger.info("Boltz detectado correctamente.")
    except FileNotFoundError:
        logger.error("Boltz no encontrado en PATH. Instalar con: pip install boltz[cuda]")
        sys.exit(1)
    
    # Verificar directorio de trabajo
    if not Config.WORKDIR.exists():
        logger.error(f"Directorio de trabajo no existe: {Config.WORKDIR}")
        logger.error("Por favor, crear el directorio y colocar los archivos:")
        logger.error(f"  - {Config.PROTEIN_PDB}")
        logger.error(f"  - {Config.LIGAND_PDB}")
        logger.error(f"  - {Config.CIF_FILE}")
        sys.exit(1)
    
    # Ejecutar benchmark
    benchmark = BoltzBenchmark()
    benchmark.run_benchmark()
    
    print("\n✓ Benchmark completado exitosamente.")


if __name__ == "__main__":
    main()
