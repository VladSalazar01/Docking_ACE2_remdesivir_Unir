@echo off
REM =============================================================================
REM INSTALADOR RAPIDO - Benchmark Boltz para Windows
REM =============================================================================
REM Este script configura el entorno conda e instala todas las dependencias
REM necesarias para ejecutar el benchmark de Boltz.
REM
REM Requisitos previos:
REM   - Miniconda/Anaconda instalado
REM   - Drivers NVIDIA actualizados
REM   - CUDA Toolkit 12.x instalado
REM
REM Uso: Ejecutar como administrador o en una terminal con permisos
REM =============================================================================

setlocal enabledelayedexpansion

echo.
echo ===============================================
echo   INSTALADOR - Benchmark Boltz
echo ===============================================
echo.

REM Verificar conda
where conda >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Conda no encontrado en PATH.
    echo Por favor instala Miniconda desde:
    echo   https://docs.conda.io/en/latest/miniconda.html
    pause
    exit /b 1
)

echo [INFO] Conda detectado correctamente.

REM Crear entorno si no existe
set ENV_NAME=boltz_env
conda info --envs | findstr /C:"%ENV_NAME%" >nul 2>&1
if %errorlevel% neq 0 (
    echo [INFO] Creando entorno conda: %ENV_NAME%
    conda create -n %ENV_NAME% python=3.10 -y
    if %errorlevel% neq 0 (
        echo [ERROR] Fallo al crear entorno conda.
        pause
        exit /b 1
    )
) else (
    echo [INFO] Entorno %ENV_NAME% ya existe.
)

REM Activar entorno
echo [INFO] Activando entorno %ENV_NAME%...
call conda activate %ENV_NAME%
if %errorlevel% neq 0 (
    echo [ERROR] No se pudo activar el entorno.
    pause
    exit /b 1
)

echo.
echo [INFO] Instalando PyTorch con soporte CUDA 12.1...
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
if %errorlevel% neq 0 (
    echo [WARN] Fallo instalacion PyTorch CUDA 12.1, intentando CUDA 12.4...
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu124
)

echo.
echo [INFO] Instalando Boltz...
pip install "boltz[cuda]" -U
if %errorlevel% neq 0 (
    echo [WARN] Intentando instalacion desde source...
    git clone https://github.com/jwohlwend/boltz.git boltz_temp
    cd boltz_temp
    pip install -e ".[cuda]"
    cd ..
)

echo.
echo [INFO] Instalando dependencias adicionales...
pip install rdkit biopython numpy pyyaml

echo.
echo [INFO] Verificando instalacion...
python -c "import torch; print(f'PyTorch: {torch.__version__}')"
python -c "import torch; print(f'CUDA disponible: {torch.cuda.is_available()}')"
python -c "import torch; print(f'GPUs: {torch.cuda.device_count()}')"
boltz --help >nul 2>&1 && echo Boltz: OK || echo Boltz: ERROR

echo.
echo ===============================================
echo   INSTALACION COMPLETADA
echo ===============================================
echo.
echo Para ejecutar el benchmark:
echo   1. Activa el entorno: conda activate %ENV_NAME%
echo   2. Navega al directorio: cd G:\Boltz_benchmark
echo   3. Ejecuta: python benchmark_boltz.py
echo.
echo Para verificar requisitos:
echo   python check_requirements.py
echo.

pause
