#!/usr/bin/env python3
"""
Test alternativo: Usar tu script original benchmark_smina_vinardo_v2.py
para reproducir el 3.187 Å
"""
import subprocess
import sys
from pathlib import Path

def test_reproduce_baseline():
    """Intentar reproducir el resultado de 3.187 Å"""
    
    base = Path(r"G:\Smina_benchmark")
    
    # Verificar que existe el script original
    original_script = base / "benchmark_smina_vinardo_v2.py"
    
    if not original_script.exists():
        print(f"✗ Original script not found: {original_script}")
        print(f"\nSearching for benchmark scripts...")
        for script in base.glob("benchmark*.py"):
            print(f"  Found: {script}")
        return
    
    print(f"✓ Found: {original_script}")
    print(f"\nRunning original benchmark to reproduce 3.187 Å...")
    print(f"This should use the EXACT configuration that gave you 3.187 Å\n")
    
    try:
        # Ejecutar el script original
        result = subprocess.run(
            [sys.executable, str(original_script)],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=str(base)
        )
        
        print(result.stdout)
        if result.stderr:
            print("STDERR:", result.stderr)
        
        print(f"\nReturn code: {result.returncode}")
        
    except subprocess.TimeoutExpired:
        print("✗ Timeout")
    except Exception as e:
        print(f"✗ Error: {e}")

if __name__ == "__main__":
    test_reproduce_baseline()