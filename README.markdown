# TFM: Análisis de Docking Molecular de ACE2 con Remdesivir

## Descripción del Proyecto

Este repositorio contiene los recursos necesarios para realizar y reproducir las simulaciones de *docking* molecular entre la proteína ACE2 (enzima convertidora de angiotensina 2) y el antiviral remdesivir, utilizando AutoDock Vina. El proyecto forma parte del Trabajo de Fin de Máster en Bioinformática y tiene como objetivo evaluar la afinidad de unión y las interacciones moleculares entre estos componentes, con el fin de explorar posibles implicaciones terapéuticas en el contexto de SARS-CoV-2.

## Estructura del Repositorio

El repositorio está organizado en las siguientes carpetas:

- **/input**: Archivos de entrada para AutoDock Vina, incluyendo estructuras de ACE2 (`ace2.pdb`, `receptor.pdbqt`), remdesivir (`remdesivir.mol2`, `ligand.pdbqt`) y configuración de la cuadrícula (`config.txt`).
- **/scripts**: Scripts en Python para preparación de receptor y ligando (`prepare_receptor.py`, `prepare_ligand.py`) y análisis de resultados (`analyze_results.py`).
- **/results**: Resultados de las simulaciones, incluyendo archivos de salida de AutoDock Vina (`docking_results.log`, `poses.pdbqt`) y visualizaciones en PyMOL (`visualization.pse`).
- **/figures**: Imágenes generadas de las interacciones moleculares y comparación de poses.
- **/docs**: Documentación adicional, como descripción de la metodología (`methods.md`) y borradores del TFM (opcional).
- **LICENSE**: Licencia MIT para el uso académico de los recursos.
- **requirements.txt**: Lista de dependencias necesarias para ejecutar las simulaciones.

## Requisitos

Para reproducir las simulaciones, se requiere instalar las siguientes herramientas en Ubuntu 22.04:

- **AutoDock Vina** (versión >= 1.2.3)
- **MGLTools/AutoDockTools** (versión 1.5.7)
- **PyMOL** (versión open-source)
- **Python** (>= 3.8) con bibliotecas: `openbabel`, `numpy`, `pandas`
- **Open Babel** (>= 3.1.1)

## Instalación

Ejecutar los siguientes comandos en Ubuntu 22.04 para configurar el entorno:

```bash
# Actualizar sistema
sudo apt update && sudo apt upgrade -y

# Instalar AutoDock Vina
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.3/vina_1.2.3_linux_x86_64 -O /usr/local/bin/vina
chmod +x /usr/local/bin/vina

# Instalar MGLTools (AutoDockTools)
wget http://mgltools.scripps.edu/downloads/tars/releases/REL1.5.7/mgltools_Linux-x86_64_1.5.7.tar.gz
tar -xvzf mgltools_Linux-x86_64_1.5.7.tar.gz
sudo mv mgltools_Linux-x86_64_1.5.7 /opt/mgltools
sudo ln -s /opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py /usr/local/bin/prepare_receptor4.py
sudo ln -s /opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py /usr/local/bin/prepare_ligand4.py
echo 'export MGLTOOLS_HOME=/opt/mgltools' >> ~/.bashrc
echo 'export PATH=$PATH:$MGLTOOLS_HOME/bin' >> ~/.bashrc
source ~/.bashrc

# Instalar PyMOL y dependencias
sudo apt install -y python3-pip python3-opengl freeglut3-dev libpng-dev libfreetype6-dev
pip3 install pymol

# Instalar Open Babel y bibliotecas de Python
sudo apt install -y openbabel
pip3 install numpy pandas openbabel
```

## Instrucciones de Uso

1. **Clonar el repositorio**:
   ```bash
   git clone https://github.com/VladSalazar01/Docking_ACE2_remdesivir_Unir
   cd tfm_docking_ace2_remdesivir
   ```

2. **Preparar receptor y ligando**:
   - Ejecutar scripts de preparación:
     ```bash
     python3 scripts/prepare_receptor.py -r input/ace2.pdb -o input/receptor.pdbqt
     python3 scripts/prepare_ligand.py -l input/remdesivir.mol2 -o input/ligand.pdbqt
     ```

3. **Ejecutar docking**:
   - Usar AutoDock Vina con el archivo de configuración:
     ```bash
     vina --config input/config.txt --receptor input/receptor.pdbqt --ligand input/ligand.pdbqt --out results/poses.pdbqt --log results/docking_results.log
     ```

4. **Analizar resultados**:
   - Procesar resultados con el script de análisis:
     ```bash
     python3 scripts/analyze_results.py results/poses.pdbqt
     ```
   - Visualizar poses en PyMOL:
     ```bash
     pymol results/visualization.pse
     ```

## Metodología

El proyecto utiliza la estructura de ACE2 obtenida del Protein Data Bank (PDB ID: 1R42) y remdesivir de la base de datos ZINC. La preparación incluye limpieza de la proteína, adición de hidrógenos y conversión a formato PDBQT. Las simulaciones de *docking* se configuran con AutoDock Vina, empleando un algoritmo genético para explorar poses de unión y calcular energías de afinidad. Los resultados se analizan mediante métricas como energía libre de unión y RMSD, con visualización de interacciones en PyMOL.

## Licencia

Este proyecto está licenciado bajo la [Licencia MIT](LICENSE), permitiendo su uso y distribución para fines académicos.

## Contacto

Para consultas relacionadas con el proyecto, contactar al autor a través del correo institucional o mediante issues en este repositorio.
