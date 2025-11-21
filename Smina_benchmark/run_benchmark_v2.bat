@echo off
REM ============================================================================
REM Smina Multi-Scoring Benchmark Launcher
REM ============================================================================

setlocal enabledelayedexpansion

set "WORK_DIR=G:\Smina_benchmark"

cd /d "%WORK_DIR%"

echo.
echo ╔═══════════════════════════════════════════════════════════════════╗
echo ║              SMINA BENCHMARK LAUNCHER                             ║
echo ╚═══════════════════════════════════════════════════════════════════╝
echo.

REM Check Python
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python not found
    pause
    exit /b 1
)

REM Check Docker
docker --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Docker not found or not running
    pause
    exit /b 1
)

echo [OK] Python and Docker detected
echo.

REM Menu
echo SELECT BENCHMARK:
echo.
echo 1. Master Suite (interactive menu)
echo 2. Multi-Scoring Comparison
echo 3. Aggressive Optimization
echo 4. Diagnostics Only
echo 5. Analysis Only
echo.
echo 0. Exit
echo.

set /p "CHOICE=Select (0-5): "

if "%CHOICE%"=="0" (
    echo Exiting...
    exit /b 0
)

if "%CHOICE%"=="1" (
    echo.
    echo Starting Master Suite...
    python master_benchmark.py
)

if "%CHOICE%"=="2" (
    echo.
    echo Starting Multi-Scoring Comparison...
    python benchmark_smina_vinardo_v2.py
)

if "%CHOICE%"=="3" (
    echo.
    echo Starting Aggressive Optimization...
    python benchmark_aggressive_optimization.py
)

if "%CHOICE%"=="4" (
    echo.
    echo Running Diagnostics...
    python diagnostic_tool.py
)

if "%CHOICE%"=="5" (
    echo.
    echo Running Analysis...
    python analyze_multiscoring.py
)

echo.
echo ════════════════════════════════════════════════════════════════════
echo COMPLETE
echo ════════════════════════════════════════════════════════════════════
echo.

pause
