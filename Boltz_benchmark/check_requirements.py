#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
=============================================================================
VERIFICADOR DE REQUISITOS - Benchmark Boltz
=============================================================================

Script para verificar que todos los requisitos están instalados correctamente
antes de ejecutar el benchmark de Boltz.

Uso:
    python check_requirements.py
"""

import sys
import os
import subprocess
from pathlib import Path

# Códigos de color para terminal
class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    END = '\033[0m'
    BOLD = '\033[1m'

def print_header():
    """Imprimir encabezado."""
    print(f"""
{Colors.BLUE}{Colors.BOLD}╔═══════════════════════════════════════════════════════════════╗
║     VERIFICADOR DE REQUISITOS - Benchmark Boltz               ║
╚═══════════════════════════════════════════════════════════════╝{Colors.END}
    """)

def check_pass(msg):
    print(f"  {Colors.GREEN}✓{Colors.END} {msg}")

def check_fail(msg, suggestion=""):
    print(f"  {Colors.RED}✗{Colors.END} {msg}")
    if suggestion:
        print(f"    {Colors.YELLOW}→ {suggestion}{Colors.END}")

def check_warn(msg):
    print(f"  {Colors.YELLOW}⚠{Colors.END} {msg}")

def check_python_version():
    """Verificar versión de Python."""
    print(f"\n{Colors.BOLD}[1/7] Verificando Python...{Colors.END}")
    
    version = sys.version_info
    if version.major == 3 and version.minor >= 9:
        check_pass(f"Python {version.major}.{version.minor}.{version.micro}")
        return True
    else:
        check_fail(
            f"Python {version.major}.{version.minor}.{version.micro} (requiere >= 3.9)",
            "Instalar Python 3.9+ o usar conda: conda create -n boltz python=3.10"
        )
        return False

def check_cuda():
    """Verificar CUDA."""
    print(f"\n{Colors.BOLD}[2/7] Verificando CUDA...{Colors.END}")
    
    # Verificar nvidia-smi
    try:
        result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
        if result.returncode == 0:
            # Extraer versión de CUDA
            for line in result.stdout.split('\n'):
                if 'CUDA Version' in line:
                    cuda_version = line.split('CUDA Version:')[1].split()[0]
                    check_pass(f"CUDA {cuda_version} disponible")
                    break
            
            # Contar GPUs
            gpu_count = result.stdout.count('RTX') + result.stdout.count('GTX') + result.stdout.count('Tesla')
            if gpu_count == 0:
                # Contar de otra forma
                gpu_count = result.stdout.count('MiB |') // 2
            
            if gpu_count >= 2:
                check_pass(f"{gpu_count} GPUs detectadas (multi-GPU habilitado)")
            elif gpu_count == 1:
                check_warn("1 GPU detectada (para mejor rendimiento usar 2 GPUs)")
            
            return True
        else:
            check_fail("nvidia-smi falló", "Verificar instalación de drivers NVIDIA")
            return False
    except FileNotFoundError:
        check_fail("nvidia-smi no encontrado", "Instalar drivers NVIDIA")
        return False

def check_pytorch():
    """Verificar PyTorch con CUDA."""
    print(f"\n{Colors.BOLD}[3/7] Verificando PyTorch...{Colors.END}")
    
    try:
        import torch
        check_pass(f"PyTorch {torch.__version__} instalado")
        
        if torch.cuda.is_available():
            check_pass("CUDA disponible en PyTorch")
            n_gpus = torch.cuda.device_count()
            for i in range(n_gpus):
                name = torch.cuda.get_device_name(i)
                mem = torch.cuda.get_device_properties(i).total_memory / (1024**3)
                check_pass(f"GPU {i}: {name} ({mem:.1f} GB)")
            return True
        else:
            check_fail(
                "CUDA NO disponible en PyTorch",
                "Reinstalar PyTorch con CUDA: pip install torch --index-url https://download.pytorch.org/whl/cu121"
            )
            return False
    except ImportError:
        check_fail(
            "PyTorch no instalado",
            "Instalar: pip install torch --index-url https://download.pytorch.org/whl/cu121"
        )
        return False

def check_boltz():
    """Verificar instalación de Boltz."""
    print(f"\n{Colors.BOLD}[4/7] Verificando Boltz...{Colors.END}")
    
    try:
        result = subprocess.run(['boltz', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            version = result.stdout.strip() or "versión desconocida"
            check_pass(f"Boltz instalado ({version})")
            return True
        else:
            # Intentar con --help
            result = subprocess.run(['boltz', '--help'], capture_output=True, text=True)
            if result.returncode == 0:
                check_pass("Boltz instalado")
                return True
            else:
                check_fail("Boltz no funciona correctamente")
                return False
    except FileNotFoundError:
        check_fail(
            "Boltz no instalado",
            "Instalar: pip install 'boltz[cuda]'"
        )
        return False

def check_rdkit():
    """Verificar RDKit."""
    print(f"\n{Colors.BOLD}[5/7] Verificando RDKit...{Colors.END}")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolAlign
        check_pass("RDKit instalado (cálculo MCS-RMSD disponible)")
        return True
    except ImportError:
        check_warn(
            "RDKit no instalado (se usará método RMSD simple)"
        )
        print(f"    {Colors.YELLOW}→ Para mejor precisión: pip install rdkit{Colors.END}")
        return True  # No es crítico

def check_biopython():
    """Verificar BioPython."""
    print(f"\n{Colors.BOLD}[6/7] Verificando BioPython...{Colors.END}")
    
    try:
        from Bio.PDB import PDBParser, MMCIFParser
        check_pass("BioPython instalado")
        return True
    except ImportError:
        check_fail(
            "BioPython no instalado",
            "Instalar: pip install biopython"
        )
        return False

def check_workdir():
    """Verificar directorio de trabajo y archivos."""
    print(f"\n{Colors.BOLD}[7/7] Verificando archivos de trabajo...{Colors.END}")
    
    workdir = Path(r"G:\Boltz_benchmark")
    
    if not workdir.exists():
        check_fail(
            f"Directorio no existe: {workdir}",
            f"Crear directorio: mkdir {workdir}"
        )
        return False
    
    check_pass(f"Directorio existe: {workdir}")
    
    # Verificar archivos
    files_to_check = [
        ("9FMM_protein.pdb", "Proteína ACE2"),
        ("ligand_A1IDX.pdb", "Ligando experimental"),
        ("9FMM.cif", "Archivo CIF (opcional)")
    ]
    
    all_present = True
    for filename, description in files_to_check:
        filepath = workdir / filename
        if filepath.exists():
            size = filepath.stat().st_size / 1024  # KB
            check_pass(f"{filename} ({size:.1f} KB) - {description}")
        else:
            if "opcional" in description.lower():
                check_warn(f"{filename} no encontrado ({description})")
            else:
                check_fail(f"{filename} no encontrado ({description})")
                all_present = False
    
    return all_present

def main():
    """Función principal."""
    print_header()
    
    results = {
        'python': check_python_version(),
        'cuda': check_cuda(),
        'pytorch': check_pytorch(),
        'boltz': check_boltz(),
        'rdkit': check_rdkit(),
        'biopython': check_biopython(),
        'workdir': check_workdir()
    }
    
    # Resumen
    print(f"\n{Colors.BOLD}{'='*60}{Colors.END}")
    print(f"{Colors.BOLD}RESUMEN{Colors.END}")
    print(f"{'='*60}")
    
    critical_ok = all([
        results['python'],
        results['cuda'],
        results['pytorch'],
        results['boltz']
    ])
    
    files_ok = results['workdir']
    
    if critical_ok and files_ok:
        print(f"\n{Colors.GREEN}{Colors.BOLD}✓ TODOS LOS REQUISITOS SATISFECHOS{Colors.END}")
        print("\nPuedes ejecutar el benchmark con:")
        print(f"  {Colors.BLUE}python benchmark_boltz.py{Colors.END}")
    elif critical_ok:
        print(f"\n{Colors.YELLOW}{Colors.BOLD}⚠ REQUISITOS DE SOFTWARE OK - FALTAN ARCHIVOS{Colors.END}")
        print("\nColoca los archivos requeridos en G:\\Boltz_benchmark\\")
    else:
        print(f"\n{Colors.RED}{Colors.BOLD}✗ FALTAN REQUISITOS CRÍTICOS{Colors.END}")
        print("\nRevisa los errores arriba e instala los componentes faltantes.")
    
    print(f"\n{'='*60}\n")
    
    return 0 if (critical_ok and files_ok) else 1

if __name__ == "__main__":
    sys.exit(main())