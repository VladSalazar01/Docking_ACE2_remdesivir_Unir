# Cribado Virtual Basado en la Estructura de ACE2

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Docker](https://img.shields.io/badge/docker-required-blue.svg)](https://www.docker.com/)

> **Trabajo de Fin de Máster** — Máster Universitario en Bioinformática  
> Universidad Internacional de La Rioja (UNIR)

## Descripción del Proyecto

Este repositorio contiene los scripts, pipelines y datos de configuración empleados en el estudio comparativo de herramientas computacionales de docking molecular de libre acceso para el cribado virtual de moduladores de la enzima convertidora de angiotensina 2 (ACE2). El trabajo aborda la validación de protocolos mediante re-docking de la estructura cristalográfica 9FMM (ACE2 humana en complejo con F-MLN-4760) y la posterior validación retrospectiva para aplicación en descubrimiento de fármacos dirigidos a enfermedades cardiovasculares.

## Objetivo General

Evaluar comparativamente el rendimiento de herramientas computacionales de docking molecular de libre acceso mediante re-docking de la estructura cristalográfica ACE2-inhibidor (PDB: 9FMM), con el propósito de establecer protocolos de cribado virtual validados para la identificación de moduladores farmacológicos de ACE2 como diana terapéutica cardiovascular.

---

## Estructura del Repositorio

```
├── Autodock_Vina_benchmark - rigido/   # Benchmark AutoDock Vina (receptor rígido)
│   ├── benchmark_vina_docker.bat       # Script principal de ejecución
│   ├── calculate_rmsd.py               # Cálculo de RMSD
│   └── resultados/                     # Outputs y logs
│
├── Gnina_benchmark/                    # Benchmark Gnina (CNN-scoring)
│   ├── benchmark_gnina_dual_site.py    # Pipeline dual-site automatizado
│   ├── calculate_rmsd_mcs.py           # RMSD con método MCS topology-aware
│   └── resultados/                     # Outputs SDF y análisis
│
├── Smina_benchmark/                    # Benchmark Smina/Vinardo
│   ├── benchmark_smina_vinardo.py      # Pipeline con múltiples funciones de scoring
│   ├── grid_kabsch.py                  # Optimización paramétrica
│   └── multiseed_test.py               # Evaluación de reproducibilidad
│
├── Dock6_benchmark_re/                 # Benchmark DOCK6 (WSL2/Ubuntu)
│   ├── run_dock6_benchmark.sh          # Script de ejecución Linux
│   ├── analyze_dock6_results.py        # Análisis de resultados
│   └── dock6.in                        # Archivo de configuración
│
├── Boltz_benchmark/                    # Benchmark Boltz-2 (predicción ab initio)
│   ├── benchmark_boltz_redocking.py    # Pipeline de predicción estructural
│   └── input.yaml                      # Configuración YAML
│
├── Cribado_VC_ACE2/                    # Validación retrospectiva y cribado virtual
│   ├── ACE1_validation/                # Validación cruzada con DUD-E ACE1
│   ├── ACE2_validation/                # Validación directa con dataset ACE2
│   ├── run_gnina_validation.py         # Ejecución de docking Docker
│   └── calculate_metrics.py            # Cálculo AUC-ROC, EF, BEDROC
│
├── ace2_remdesivir_docking.sh          # Script legacy (exploración inicial)
└── README.md                           # Este archivo
```

---

## Resultados del Benchmarking

### Comparación de Precisión Geométrica (Re-docking 9FMM)

| Programa | RMSD (Å) | Desv. Est. | N runs | Success ≤2Å | Tiempo/run | Aceleración |
|----------|:--------:|:----------:|:------:|:-----------:|:----------:|:-----------:|
| **DOCK6** | **0.473** | 0.009 | 5 | **100%** | ~22.7 min | CPU |
| **Gnina** | **1.544** | 0.363 | 17 | **88.2%** | ~25 s | GPU |
| Smina/Vinardo | 3.320 | 0.090 | 20 | 0% | ~155 s | CPU |
| AutoDock Vina | 4.463† | 0.086 | 5 | 0% | ~45.7 s | CPU |
| Boltz-2 | 5.240 | 1.340 | 10 | 0% | ~96 s | GPU |

> † RMSD de AutoDock Vina afectado por configuración de caja de docking subóptima.  
> Criterio de éxito: RMSD ≤ 2.0 Å respecto a la conformación cristalográfica.

### Validación Retrospectiva (Capacidad Discriminativa)

| Diana | Métrica | Gnina | Umbral | Evaluación |
|-------|---------|:-----:|:------:|:----------:|
| ACE1 (DUD-E) | AUC-ROC | 0.641 | ≥0.70 | Moderado |
| ACE1 (DUD-E) | EF 1% | 6.82 | ≥10 | Aceptable |
| **ACE2** | **AUC-ROC** | **0.748** | ≥0.70 | **Bueno** |
| **ACE2** | **EF 1%** | **5.26** | -- | Aceptable |
| ACE2 | BEDROC (α=20) | 0.298 | ≥0.30 | Moderado |

---

## Requisitos del Sistema

### Hardware Empleado
- **CPU:** AMD Ryzen Threadripper 3690X (24 cores / 48 threads)
- **GPU:** 2× NVIDIA RTX 3060 (12 GB VRAM cada una)
- **RAM:** 128 GB DDR4
- **Almacenamiento:** SSD NVMe (recomendado para I/O intensivo)

### Software
- **Sistema Operativo:** Windows 10 x64 / WSL2 Ubuntu 22.04
- **Docker Desktop:** Para contenedores de Gnina y otros programas
- **Python:** 3.9+
- **Dependencias principales:**
  - NumPy ≥1.23
  - RDKit ≥2022.09
  - SciPy ≥1.9
  - Matplotlib ≥3.6

### Imágenes Docker Utilizadas
```bash
# Gnina (CNN-scoring docking)
docker pull gnina/gnina:latest

# AutoDock Vina
docker pull ccsb/autodock-vina:latest
```

---

## Instalación y Configuración

### 1. Clonar el Repositorio
```bash
git clone https://github.com/[usuario]/ace2-virtual-screening.git
cd ace2-virtual-screening
```

### 2. Crear Entorno Virtual (Python)
```bash
python -m venv venv
source venv/bin/activate  # Linux/macOS
# o
.\venv\Scripts\activate   # Windows

pip install -r requirements.txt
```

### 3. Verificar Docker
```bash
docker --version
docker run --rm gnina/gnina:latest gnina --version
```

### 4. Configurar DOCK6 (WSL2)
```bash
# Desde Windows PowerShell
wsl --install -d Ubuntu-22.04

# Dentro de WSL
sudo apt update && sudo apt install gfortran gcc make
cd /mnt/g/Dock6_benchmark_re
tar -xzf dock6-latest.tgz
```

---

## Uso Básico

### Benchmark de Gnina (Recomendado)
```bash
cd Gnina_benchmark
python benchmark_gnina_dual_site.py --runs 17 --exhaustiveness 64
```

### Benchmark de DOCK6 (WSL2)
```bash
cd /mnt/g/Dock6_benchmark_re
bash run_dock6_benchmark.sh
python3 analyze_dock6_results.py
```

### Validación Retrospectiva ACE2
```bash
cd Cribado_VC_ACE2
python run_gnina_validation.py --receptor 9FMM_prepared.pdbqt --actives activos_ace2.sdf --decoys decoys_ace2.sdf
python calculate_metrics.py --results docking_results.csv
```

---

## Metodología

### Preparación de Estructuras
1. **Receptor (9FMM):** Eliminación de aguas, adición de hidrógenos polares, asignación de cargas Gasteiger
2. **Ligando (F-MLN-4760):** Extracción del PDB, minimización con MMFF94, generación de conformaciones

### Cálculo de RMSD
Se implementó el método **MCS (Maximum Common Substructure)** para cálculos topology-aware que evitan problemas de orden de átomos entre formatos:

```python
from rdkit.Chem import AllChem, rdMolAlign

rmsd = rdMolAlign.GetBestRMS(mol_ref, mol_pred)
```

### Métricas de Validación
- **AUC-ROC:** Área bajo la curva ROC para discriminación activos/decoys
- **EF (Enrichment Factor):** Factor de enriquecimiento al 1%, 5% y 10%
- **BEDROC:** Boltzmann-Enhanced Discrimination of ROC (α=20)

---

## Estructura Cristalográfica de Referencia

| Parámetro | Valor |
|-----------|-------|
| **PDB ID** | 9FMM |
| **Resolución** | 2.15 Å |
| **Proteína** | ACE2 humana (ectodominio) |
| **Ligando** | F-MLN-4760 (inhibidor) |
| **Centro del Grid** | (42.0, 7.0, 23.0) Å |
| **Tamaño del Grid** | 26 × 26 × 26 Å |

---

## Publicaciones y Referencias

### Programas de Docking
1. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking. *J Comput Chem.* 2010;31(2):455-461.
2. McNutt AT, et al. GNINA 1.0: molecular docking with deep learning. *J Cheminform.* 2021;13(1):43.
3. Allen WJ, et al. DOCK 6: Impact of new features and current docking performance. *J Comput Chem.* 2015;36(15):1132-1156.
4. Wohlwend J, et al. Boltz-1: Democratizing Biomolecular Interaction Modeling. 2024. bioRxiv.

### Datasets de Validación
5. Mysinger MM, et al. Directory of Useful Decoys, Enhanced (DUD-E). *J Med Chem.* 2012;55(14):6582-6594.

---

## Licencia

Este proyecto está licenciado bajo la Licencia MIT. Consulte el archivo [LICENSE](LICENSE) para más detalles.

---

## Autor

**Vladimir Oswaldo Salazar Córdova**  
Máster Universitario en Bioinformática  
Universidad Internacional de La Rioja (UNIR)

**Director:** Alberto Robles Loaiza

---

## Agradecimientos

- NVIDIA por soporte de GPU computing
- Comunidad de RDKit y OpenBabel

---

*Última actualización: Febrero 2026*
