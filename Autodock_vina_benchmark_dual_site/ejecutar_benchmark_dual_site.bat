@echo off
chcp 65001 >nul 2>&1
title AutoDock Vina - Benchmark Dual-Site (Corrección Retrospectiva)

echo ============================================================
echo   AUTODOCK VINA - BENCHMARK DUAL-SITE
echo   Corrección retrospectiva del centro de caja
echo ============================================================
echo.
echo   Centro INCORRECTO anterior (MID): (12.47, 15.46, 49.91)
echo   Centro Fragment 1 (corregido):    (16.48, 15.24, 25.54)
echo   Centro Fragment 2 (corregido):    (8.44,  15.48, 74.27)
echo.

:: Verificar Python
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python no encontrado. Activar conda environment.
    pause
    exit /b 1
)

:: Directorio de trabajo
cd /d "G:\Autodock_vina_benchmark"

:: Número de runs
set /p N_RUNS="Número de runs por sitio [5]: "
if "%N_RUNS%"=="" set N_RUNS=5

echo.
echo Ejecutando %N_RUNS% runs por sitio (%N_RUNS% x 2 = total)...
echo.

python benchmark_vina_dual_site.py %N_RUNS%

echo.
echo ============================================================
echo   Benchmark completado.
echo   Revisa la carpeta resultados_dual_site_* para los logs.
echo ============================================================
pause
