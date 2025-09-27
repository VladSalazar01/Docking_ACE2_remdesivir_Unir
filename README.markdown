# TFM: Análisis Comparativo de Docking Molecular de ACE2 con Remdesivir

![Docking Visualization](figures/remdesivir_ace2_pose.png)

## 📖 Descripción del Proyecto

Este repositorio contiene el pipeline bioinformático completo para el **Trabajo de Fin de Máster (TFM)** en Bioinformática, titulado **"Análisis Comparativo de Docking Molecular entre la Proteína ACE2 y Remdesivir"**. El proyecto evalúa la afinidad de unión e interacciones moleculares del antiviral remdesivir con la enzima convertidora de angiotensina 2 (ACE2), una diana clave en SARS-CoV-2, utilizando **AutoDock Vina** como herramienta principal. Para enriquecer el análisis, se implementa un enfoque de **screening virtual (VS)** validado, que incluye:

- Búsqueda de ligandos bioactivos (IC50 ≤10 μM o <100 nM) en bases de datos como ZINC, DrugBank y ChEMBL, usando **RDKit** fingerprints (Morgan, ECFP, MACCS).
- Generación de decoys con **DUDE-Z** para benchmarking.
- Refinamiento cuántico semi-empírico con **MOPAC** (PM7, MOZYME para complejos).
- Comparación de herramientas de docking open-source: **Vina, Smina, GNINA, VINARDO, AutoDock4 (AD4), DOCK6, BOLTZ2**.
- Análisis estadístico (ROC-AUC, correlaciones de scores) y visualización en **PyMOL**.

El objetivo es explorar implicaciones terapéuticas de remdesivir como posible modulador de ACE2, validando la fiabilidad de las herramientas de docking en un contexto de VS. Este pipeline híbrido (ligando-based + estructura-based) aporta rigor científico, reproducibilidad y potencial para extensiones (e.g., simulaciones de dinámica molecular). Los resultados incluyen métricas de afinidad (kcal/mol), RMSD, y visualizaciones de interacciones (H-bonds, hidrofóbicas), con un enfoque comparativo que eleva el TFM a un estándar publicable.

## 🎯 Objetivos

- Evaluar la interacción molecular ACE2-remdesivir en el sitio de unión RBD (receptor-binding domain).
- Comparar el desempeño de múltiples herramientas de docking para robustez.
- Validar el pipeline de VS mediante decoys y métricas estadísticas (e.g., enrichment factor >10).
- Proporcionar un flujo reproducible documentado en Markdown y Jupyter Notebooks.

## 🚀 Características

- **Pipeline Completo**: Desde búsqueda de ligandos hasta análisis comparativo y visualización.
- **Validación Robusta**: Uso de decoys y métricas ROC-AUC para evaluar especificidad.
- **Cuantificación Cuántica**: Refinamiento con MOPAC para geometrías y cargas precisas.
- **Comparación Multi-Herramienta**: Evaluación de Vina, Smina, GNINA, etc., con correlaciones de scores.
- **Reproducibilidad**: Entorno Conda, scripts modulares, y documentación detallada.
- **Visualización Avanzada**: Poses e interacciones en PyMOL; gráficos de resultados en Matplotlib/Seaborn.

## 📂 Estructura del Repositorio

```plaintext
tfm_docking_ace2_remdesivir/
├── input/                   # Archivos de entrada (ACE2: 6M0J.pdb, remdesivir.smi, config_vina.txt)
├── ligands/                 # Ligandos bioactivos (SMILES, PDB; ~10 con IC50 ≤10 μM)
├── decoys/                  # Decoys generados por DUDE-Z (~50 por ligando)
├── mopac/                   # Inputs/outputs de refinamiento QM (e.g., remdesivir.mop, charges.out)
├── scripts/                 # Scripts Python para el pipeline
│   ├── prepare_ligands.py   # Búsqueda y preparación de ligandos con RDKit/OpenBabel
│   ├── generate_decoys.py   # Generación de decoys desde DUDE-Z
│   ├── refine_mopac.py      # Refinamiento QM con MOPAC (PM7/MOZYME)
│   ├── run_dockings.py      # Docking comparativo (Vina, Smina, etc.)
│   ├── analyze_results.py    # Análisis estadístico (ROC, RMSD, scores)
│   ├── prepare_receptor.py   # Preparación de ACE2 para docking
│   └── prepare_ligand.py     # Conversión de ligandos a PDBQT
├── dockings/                # Resultados por herramienta (/vina, /smina, etc.; poses.pdbqt, logs)
├── results/                 # Outputs consolidados (affinity_scores.csv, roc_auc.png)
├── figures/                 # Visualizaciones (PyMOL poses, plots Matplotlib)
├── docs/                    # Documentación (methods.md, tfm_draft.md)
├── analysis/                # Notebooks Jupyter (e.g., rdkit_analysis.ipynb)
├── environment.yml          # Entorno Conda para reproducibilidad
├── requirements.txt         # Dependencias Python
└── LICENSE                  # Licencia MIT
```

## 🛠 Requisitos

- **Sistema Operativo**: Ubuntu 24.04+ (compatible con 2025) o WSL2 en Windows.
- **Hardware**: CPU multi-core (≥8 cores); GPU NVIDIA (recomendado para GNINA/BOLTZ2).
- **Dependencias**:
  - Python (≥3.10): RDKit, NumPy, Pandas, Matplotlib, Seaborn, OpenBabel (wrapper).
  - Herramientas de docking: AutoDock Vina (≥1.2.7), Smina, GNINA, VINARDO, AutoDock4, DOCK6, BOLTZ2.
  - MGLTools/AutoDockTools (1.5.7).
  - MOPAC (≥2023 para MOZYME).
  - PyMOL (open-source ≥3.1.6).
  - Open Babel (≥3.1.1).
- **Acceso Web**: ZINC, ChEMBL, DrugBank (ligandos); DUDE-Z (decoys).

## 🔧 Instalación

Se recomienda usar **Conda** para un entorno reproducible. Sigue estos pasos en Ubuntu 24.04+:

```bash
# 1. Clonar el repositorio
git clone https://github.com/VladSalazar01/Docking_ACE2_remdesivir_Unir
cd tfm_docking_ace2_remdesivir

# 2. Instalar Miniconda (si no está instalado)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b
source ~/.bashrc

# 3. Crear y activar entorno Conda
conda env create -f environment.yml
conda activate tfm_docking

# 4. Instalar AutoDock Vina
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_linux_x86_64 -O ~/bin/vina
chmod +x ~/bin/vina
export PATH=$PATH:~/bin

# 5. Instalar MGLTools
wget http://mgltools.scripps.edu/downloads/tars/releases/REL1.5.7/mgltools_Linux-x86_64_1.5.7.tar.gz
tar -xvzf mgltools_Linux-x86_64_1.5.7.tar.gz
sudo mv mgltools_Linux-x86_64_1.5.7 /opt/mgltools
echo 'export MGLTOOLS_HOME=/opt/mgltools' >> ~/.bashrc
echo 'export PATH=$PATH:$MGLTOOLS_HOME/bin:/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24' >> ~/.bashrc
source ~/.bashrc

# 6. Instalar PyMOL y Open Babel
pip install pymol-open-source-whl
sudo apt update && sudo apt install -y openbabel

# 7. Herramientas adicionales (Smina, GNINA, etc.)
# Smina: conda install -c conda-forge smina
# GNINA: git clone https://github.com/gnina/gnina.git; cd gnina; mkdir build; cd build; cmake ..; make -j4
# VINARDO, AD4, DOCK6, BOLTZ2: Ver docs/INSTALL_TOOLS.md
```

**Nota**: MOPAC y RDKit están en `environment.yml`. Para Google Colab, usa `notebooks/labodock_colab.ipynb` (limitado para MOPAC/GNINA).

## 🚀 Instrucciones de Uso

Ejecuta el pipeline completo con `bash run_pipeline.sh` o paso a paso:

1. **Preparar receptor** (ACE2, PDB ID: 6M0J):
   ```bash
   python scripts/prepare_receptor.py -r input/ace2.pdb -o input/receptor.pdbqt
   ```

2. **Buscar y preparar ligandos** (remdesivir + ~10 similares):
   ```bash
   python scripts/prepare_ligands.py --query "ACE2 inhibitors IC50 <10uM" --output ligands/
   # Usa RDKit para fingerprints; OpenBabel para SMILES → PDB
   ```

3. **Generar decoys** (DUDE-Z):
   ```bash
   # Subir SMILES manualmente a dudez.docking.org; luego:
   python scripts/generate_decoys.py --smiles ligands/remdesivir.smi --decoys_dir decoys/
   ```

4. **Refinar con MOPAC** (geometrías y cargas):
   ```bash
   python scripts/refine_mopac.py --input ligands/*.pdb --method PM7 --mozyme True --output mopac/
   ```

5. **Ejecutar dockings comparativos**:
   ```bash
   python scripts/run_dockings.py --config input/config.txt --tools all --output dockings/
   # Grid box: ~20x20x20 Å en sitio RBD; exhaustiveness=8-32
   ```

6. **Analizar resultados**:
   ```bash
   python scripts/analyze_results.py --dockings dockings/ --decoys decoys/ --output results/
   # Genera: affinity_scores.csv, roc_auc.png, RMSD stats
   ```

7. **Visualizar**:
   ```bash
   pymol scripts/visualize.py  # Carga poses e interacciones
   ```

**Output Ejemplo**:
- `results/affinity_scores.csv`: `[Tool, Ligand, Binding_Energy_kcal/mol, RMSD]`
- `figures/roc_auc.png`: Curva ROC para validación de VS.

## 📝 Metodología

El pipeline sigue un flujo bioinformático robusto:

1. **Estructuras**: ACE2 (PDB 6M0J, resolución 2.5 Å, limpio con PyMOL). Remdesivir y ligandos similares (ChEMBL/ZINC, filtrados por IC50 y Tanimoto >0.7).
2. **Ligandos**: SMILES → PDB con OpenBabel; fingerprints RDKit (Morgan/ECFP).
3. **Decoys**: ~50 por ligando desde DUDE-Z (MW/logP matched).
4. **Refinamiento QM**: MOPAC PM7 (geometrías); MOZYME para cargas de complejos.
5. **Docking**: Multi-herramienta (Vina como base); métricas: ΔG, RMSD <2 Å.
6. **Análisis**: ROC-AUC (>0.8), correlaciones Pearson de scores (>0.8), visualización de H-bonds/hidrofóbicas en PyMOL.

**Contribuciones al TFM**:
- Comparación multi-herramienta para robustez.
- Validación estadística de VS.
- Análisis cuántico para precisión electrostática.
- Documentación reproducible para publicación académica.

## 📊 Resultados Esperados

- **Afinidad de Remdesivir**: ΔG ~-7 a -9 kcal/mol (Vina); comparación con inhibidores conocidos.
- **Validación**: ROC-AUC >0.8; enrichment factor >10 para actives vs. decoys.
- **Interacciones**: H-bonds con residuos clave (e.g., Asp350, His34 en ACE2).
- **Gráficos**: Scatter plots (scores Vina vs. GNINA), visualizaciones 3D.

## 📜 Licencia

Este proyecto está bajo la **[Licencia MIT](LICENSE)**, permitiendo uso, modificación y distribución para fines académicos y de investigación no comercial.

## 🙌 Contribuciones

¡Bienvenidas las contribuciones! Por favor, revisa [CONTRIBUTING.md](docs/CONTRIBUTING.md) para directrices. Abre un issue o pull request para mejoras (e.g., scripts adicionales, optimización GPU).

## 📧 Contacto

- **Autor**: Vlad Salazar ([vosalazar26@outlook.com](mailto:vosalazar26@outlook.com))
- **Issues**: Reporta bugs o sugerencias en [GitHub Issues](https://github.com/VladSalazar01/Docking_ACE2_remdesivir_Unir/issues)
- **Colaboración**: Contacta para extensiones (e.g., MD con GROMACS) o publicaciones.

## 🌟 Agradecimientos

- Equipo de AutoDock Vina, RDKit, MOPAC, y DUDE-Z por herramientas open-source.
- Comunidad de bioinformática por recursos como ChEMBL y ZINC.
- UNIR por soporte académico.
