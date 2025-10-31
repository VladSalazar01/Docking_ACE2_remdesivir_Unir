@echo off
setlocal enabledelayedexpansion

:: Activar environment de conda
call conda activate pymol_env

:: Configuración
set "WORKDIR=G:\Autodock_vina_benchmark"
cd /d %WORKDIR% 2>nul
if errorlevel 1 (
    echo Error: No se pudo cambiar a %WORKDIR%. Verifica la unidad G:.
    pause
    exit /b
)
echo Directorio de trabajo establecido en %WORKDIR%

set "IMAGE=ghcr.io/metaphorme/vina:v1.2.5"
set "INPUT_FILE=%WORKDIR%\9FMM.cif"
set "RECEPTOR=%WORKDIR%\receptor.pdbqt"
set "LIGAND=%WORKDIR%\ligand.pdbqt"
set "CONFIG=%WORKDIR%\config.txt"
set "OUTPUT=%WORKDIR%\output.pdbqt"
set "LIG_EXP=%WORKDIR%\ligand_exp.pdb"
set "RUNS=5"
set "CPU=48"

:: Ruta a Python (conda environment) y Open Babel
set "PYTHON=python"
set "OBABEL=C:\Program Files\OpenBabel-3.1.1\obabel.exe"

:: Verificar Open Babel
if not exist "%OBABEL%" (
    echo Error: Open Babel no encontrado en %OBABEL%
    pause
    exit /b
)

:: Verificar si el archivo de entrada existe y tiene datos
if not exist "%INPUT_FILE%" (
    echo Error: %INPUT_FILE% no encontrado.
    pause
    exit /b
)
findstr "_atom_site" "%INPUT_FILE%" >nul
if errorlevel 1 (
    echo Error: %INPUT_FILE% parece corrupto o sin datos atómicos.
    echo Descarga nuevamente desde https://files.rcsb.org/download/9FMM.cif
    pause
    exit /b
)

:: Crear config.txt si no existe
if not exist "%CONFIG%" (
    echo receptor = /data/receptor.pdbqt > "%CONFIG%"
    echo ligand = /data/ligand.pdbqt >> "%CONFIG%"
    echo out = /data/output.pdbqt >> "%CONFIG%"
    echo center_x = 12.470 >> "%CONFIG%"
    echo center_y = 15.461 >> "%CONFIG%"
    echo center_z = 49.907 >> "%CONFIG%"
    echo size_x = 20 >> "%CONFIG%"
    echo size_y = 20 >> "%CONFIG%"
    echo size_z = 20 >> "%CONFIG%"
    echo exhaustiveness = 32 >> "%CONFIG%"
    echo num_modes = 1 >> "%CONFIG%"
    echo cpu = %CPU% >> "%CONFIG%"
    echo Config.txt creado en %CONFIG%.
)

:: Verificar que config.txt se creó correctamente
if not exist "%CONFIG%" (
    echo Error: No se pudo crear %CONFIG%. Verifica permisos.
    pause
    exit /b
)

:: Preparar receptor si no existe
if not exist "%RECEPTOR%" (
    echo Preparando receptor desde 9FMM_protein.pdb...
    
    if not exist "%WORKDIR%\9FMM_protein.pdb" (
        echo Error: 9FMM_protein.pdb no encontrado. Generalo primero con PyMOL.
        pause
        exit /b
    )
    
    "%OBABEL%" -i pdb "%WORKDIR%\9FMM_protein.pdb" -o pdbqt -O "%RECEPTOR%" -xr --partialcharge gasteiger
    
    if not exist "%RECEPTOR%" (
        echo Error: No se pudo generar %RECEPTOR%.
        pause
        exit /b
    )
    echo Receptor %RECEPTOR% creado exitosamente.
)

:: Preparar ligando usando Meeko si no existe
if not exist "%LIGAND%" (
    echo Preparando ligando con Meeko desde ligand_A1IDX.pdb...
    
    if not exist "%WORKDIR%\ligand_A1IDX.pdb" (
        echo Error: ligand_A1IDX.pdb no encontrado en %WORKDIR%
        pause
        exit /b
    )
    
    :: Crear script de Python para preparar ligando
    echo from meeko import MoleculePreparation, PDBQTWriterLegacy > "%WORKDIR%\prepare_ligand_temp.py"
    echo from rdkit import Chem >> "%WORKDIR%\prepare_ligand_temp.py"
    echo. >> "%WORKDIR%\prepare_ligand_temp.py"
    echo mol = Chem.MolFromPDBFile('ligand_A1IDX.pdb', removeHs=False^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo if mol is None: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     print("Error: No se pudo leer ligand_A1IDX.pdb"^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo else: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     mol = Chem.AddHs(mol, addCoords=True^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     if len(frags^) ^> 1: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo         print(f"Encontrados {len(frags^)} fragmentos. Seleccionando el mayor..."^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo         mol = max(frags, key=lambda m: m.GetNumAtoms(^)^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     preparator = MoleculePreparation(^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     mol_setups = preparator.prepare(mol^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo     for setup in mol_setups: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo         pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo         if is_ok: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo             with open('ligand.pdbqt', 'w'^) as f: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo                 f.write(pdbqt_string^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo             print("Ligando preparado: ligand.pdbqt"^) >> "%WORKDIR%\prepare_ligand_temp.py"
    echo         else: >> "%WORKDIR%\prepare_ligand_temp.py"
    echo             print(f"Error: {error_msg}"^) >> "%WORKDIR%\prepare_ligand_temp.py"
    
    %PYTHON% "%WORKDIR%\prepare_ligand_temp.py"
    del "%WORKDIR%\prepare_ligand_temp.py"
    
    if not exist "%LIGAND%" (
        echo Error: No se pudo generar %LIGAND% con Meeko.
        pause
        exit /b
    )
    
    echo Ligando %LIGAND% creado. Verificando contenido...
    findstr "ROOT" "%LIGAND%" >nul
    if errorlevel 1 (
        echo Error: %LIGAND% no contiene ROOT. Formato invalido.
        pause
        exit /b
    ) else (
        echo Ligando contiene ROOT. Formato correcto.
    )
)

:: Generar pose experimental si no existe
if not exist "%LIG_EXP%" (
    echo Generando pose experimental desde ligand_A1IDX.pdb...
    copy "%WORKDIR%\ligand_A1IDX.pdb" "%LIG_EXP%" >nul
    if not exist "%LIG_EXP%" (
        echo Error: No se pudo copiar pose experimental.
        pause
        exit /b
    )
)

:: Bucle de benchmarking
echo.
echo ========================================
echo Iniciando benchmarking con %RUNS% runs
echo ========================================
echo.

for /l %%i in (1, 1, %RUNS%) do (
    echo [Run %%i/%RUNS%] Ejecutando docking...
    
    docker run -v G:\Autodock_vina_benchmark:/data -it --rm "%IMAGE%" vina --config /data/config.txt --out /data/output_%%i.pdbqt > "%WORKDIR%\vina_log_%%i.txt" 2>&1
    
    if errorlevel 1 (
        echo [Run %%i] ERROR en docking. Revisa %WORKDIR%\vina_log_%%i.txt
    ) else (
        echo [Run %%i] Docking completado exitosamente.
        
        :: Convertir PDBQT a PDB con Open Babel (sin agregar hidrógenos)
        "%OBABEL%" -i pdbqt "%WORKDIR%\output_%%i.pdbqt" -o pdb -O "%WORKDIR%\output_%%i.pdb" -d >nul 2>&1
        
        :: Calcular RMSD usando Python + RDKit con alineamiento Kabsch
        echo from rdkit import Chem > "%WORKDIR%\calc_rmsd_temp.py"
        echo import numpy as np >> "%WORKDIR%\calc_rmsd_temp.py"
        echo. >> "%WORKDIR%\calc_rmsd_temp.py"
        echo try: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     ref = Chem.MolFromPDBFile('ligand_exp.pdb', removeHs=False^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     probe = Chem.MolFromPDBFile('output_%%i.pdb', removeHs=False^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     if ref is None: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         print("Error: No se pudo leer ligand_exp.pdb"^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     elif probe is None: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         print("Error: No se pudo leer output_%%i.pdb"^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     elif ref.GetNumAtoms(^) != probe.GetNumAtoms(^): >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         print(f"Error: Ref={ref.GetNumAtoms(^)} vs Probe={probe.GetNumAtoms(^)}"^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     else: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         ref_conf = ref.GetConformer(^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         probe_conf = probe.GetConformer(^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         n = ref.GetNumAtoms(^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         ref_coords = np.array([ref_conf.GetAtomPosition(i^) for i in range(n^)]^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         probe_coords = np.array([probe_conf.GetAtomPosition(i^) for i in range(n^)]^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         ref_center = ref_coords - ref_coords.mean(axis=0^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         probe_center = probe_coords - probe_coords.mean(axis=0^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         H = probe_center.T @ ref_center >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         U, S, Vt = np.linalg.svd(H^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         R = Vt.T @ U.T >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         if np.linalg.det(R^) ^< 0: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo             Vt[-1, :] *= -1 >> "%WORKDIR%\calc_rmsd_temp.py"
        echo             R = Vt.T @ U.T >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         probe_aligned = (R @ probe_center.T^).T >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         diff = ref_center - probe_aligned >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1^)^)^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         result = f"RMSD: {rmsd:.3f} A" >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         print(result^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo         with open('rmsd_%%i.txt', 'w'^) as f: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo             f.write(result + "\n"^) >> "%WORKDIR%\calc_rmsd_temp.py"
        echo except Exception as e: >> "%WORKDIR%\calc_rmsd_temp.py"
        echo     print(f"Error: {e}"^) >> "%WORKDIR%\calc_rmsd_temp.py"
        
        %PYTHON% "%WORKDIR%\calc_rmsd_temp.py"
        del "%WORKDIR%\calc_rmsd_temp.py"
        
        if exist "%WORKDIR%\rmsd_%%i.txt" (
            type "%WORKDIR%\rmsd_%%i.txt"
        )
    )
    echo.
)

echo ========================================
echo Benchmarking completado
echo ========================================
echo Revisa logs (vina_log_*.txt) y resultados (rmsd_*.txt) en %WORKDIR%
echo.

:: Resumen estadístico
echo Generando resumen estadistico...
echo import glob > "%WORKDIR%\summary_temp.py"
echo rmsds = [] >> "%WORKDIR%\summary_temp.py"
echo for file in sorted(glob.glob('rmsd_*.txt')): >> "%WORKDIR%\summary_temp.py"
echo     with open(file) as f: >> "%WORKDIR%\summary_temp.py"
echo         line = f.read().strip() >> "%WORKDIR%\summary_temp.py"
echo         if 'RMSD:' in line: >> "%WORKDIR%\summary_temp.py"
echo             rmsd = float(line.split('RMSD:')[1].split('A')[0].strip()) >> "%WORKDIR%\summary_temp.py"
echo             rmsds.append(rmsd) >> "%WORKDIR%\summary_temp.py"
echo if rmsds: >> "%WORKDIR%\summary_temp.py"
echo     avg = sum(rmsds) / len(rmsds) >> "%WORKDIR%\summary_temp.py"
echo     print(f"RMSD promedio: {avg:.3f} A") >> "%WORKDIR%\summary_temp.py"
echo     print(f"RMSD minimo:   {min(rmsds):.3f} A") >> "%WORKDIR%\summary_temp.py"
echo     print(f"RMSD maximo:   {max(rmsds):.3f} A") >> "%WORKDIR%\summary_temp.py"

%PYTHON% "%WORKDIR%\summary_temp.py"
del "%WORKDIR%\summary_temp.py"

pause