# ACE2 Virtual Screening Pipelines

![Python](https://img.shields.io/badge/python-3.9+-blue.svg)
![Docker](https://img.shields.io/badge/docker-required-2496ED.svg)
![Status](https://img.shields.io/badge/status-active-success.svg)
![License](https://img.shields.io/badge/license-MIT-green.svg)

Pipelines automatizados para benchmarking de docking molecular y cribado virtual sobre ACE2 (PDB: 9FMM).

## Quick Start

```bash
# Clonar repositorio
git clone https://github.com/[usuario]/ace2-virtual-screening.git
cd ace2-virtual-screening

# Instalar dependencias
pip install numpy rdkit scipy pandas matplotlib meeko

# Ejecutar benchmark Gnina (requiere Docker)
cd Gnina_benchmark
python benchmark_gnina_dual_site.py 10
```

## Estructura

```
├── Autodock_Vina_benchmark - rigido/   # AutoDock Vina dual-site
├── Gnina_benchmark/                    # Gnina con CNN-scoring (GPU)
├── Smina_benchmark/                    # Smina/Vinardo optimización
├── Dock6_benchmark_re/                 # DOCK6 (WSL2)
├── Boltz_benchmark/                    # Boltz-2 predicción ab initio
└── Cribado_VC_ACE2/                    # Validación retrospectiva
```

## Pipelines

### Gnina (GPU)

```bash
cd Gnina_benchmark
python benchmark_gnina_dual_site.py 17    # Dual-site, 17 runs
python calculate_rmsd_mcs.py ref.pdb out.pdb  # RMSD topology-aware
```

### AutoDock Vina

```bash
cd "Autodock_Vina_benchmark - rigido"
python benchmark_vina_dual_site.py --runs 5 --exhaustiveness 64
```

### Smina/Vinardo

```bash
cd Smina_benchmark
python grid_kabsch.py          # Grid search paramétrico
python multiseed_test.py       # Test reproducibilidad (20 seeds)
```

### DOCK6 (WSL2)

```bash
cd /mnt/g/Dock6_benchmark_re
bash run_dock6_benchmark.sh
python3 analyze_dock6_results.py
```

### Boltz-2

```bash
cd Boltz_benchmark
python benchmark_boltz_redocking.py gpu --runs 10
```

### Validación Retrospectiva

```bash
cd Cribado_VC_ACE2/ACE2_validation
python run_gnina_validation.py
python calculate_metrics.py    # AUC-ROC, EF, BEDROC
```

## Requisitos

- **OS:** Windows 10 + WSL2 (Ubuntu 22.04)
- **Python:** 3.9+
- **Docker:** gnina/gnina, ccsb/autodock-vina
- **Hardware:** GPU con ≥12GB VRAM (para Gnina y Boltz)

## Tecnologías

![RDKit](https://img.shields.io/badge/RDKit-2022.09-orange)
![NumPy](https://img.shields.io/badge/NumPy-1.23-blue)
![Docker](https://img.shields.io/badge/Docker-Desktop-2496ED)
![PyMOL](https://img.shields.io/badge/PyMOL-2.5-red)

## Autor

**Vladimir Oswaldo Salazar Córdova**  
MSc Bioinformática — UNIR

## Licencia

[MIT](LICENSE)
