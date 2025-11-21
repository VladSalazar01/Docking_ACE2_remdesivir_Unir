#!/usr/bin/env python3
"""
VERIFICACIÓN DE ENTORNO
=======================
Verifica que todo esté listo antes de ejecutar el pipeline
"""

import sys
import os
import subprocess
from pathlib import Path
import shutil

def check_python():
    """Verifica versión de Python"""
    print("\n1. Verificando Python...")
    version = sys.version_info
    print(f"   ✓ Python {version.major}.{version.minor}.{version.micro}")
    if version.major < 3 or (version.major == 3 and version.minor < 7):
        print("   ⚠ Recomendado: Python 3.7+")
        return False
    return True

def check_smina():
    """Verifica que Smina esté instalado y accesible"""
    print("\n2. Verificando Smina...")
    
    # Intentar ejecutar smina
    try:
        result = subprocess.run(
            ["smina", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0 or "smina" in result.stderr.lower():
            print("   ✓ Smina encontrado en PATH")
            return True
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"   Error inesperado: {e}")
    
    # Buscar smina.exe en ubicaciones comunes de Windows
    print("   ✗ Smina no está en PATH")
    print("\n   Buscando smina.exe en ubicaciones comunes...")
    
    common_locations = [
        r"C:\Program Files\Smina",
        r"C:\Program Files (x86)\Smina",
        r"C:\Smina",
        r"C:\Tools\Smina",
        os.path.expanduser("~\\AppData\\Local\\Smina"),
    ]
    
    for location in common_locations:
        smina_path = Path(location) / "smina.exe"
        if smina_path.exists():
            print(f"   ✓ Encontrado: {smina_path}")
            print(f"\n   💡 Para agregarlo al PATH:")
            print(f"      1. Copiar esta ruta: {location}")
            print(f"      2. Variables de entorno → PATH → Agregar")
            return False
    
    print("   ✗ Smina.exe no encontrado en ubicaciones comunes")
    print("\n   📥 INSTALACIÓN DE SMINA:")
    print("   1. Descargar de: https://sourceforge.net/projects/smina/")
    print("   2. Extraer smina.exe a C:\\Smina\\")
    print("   3. Agregar C:\\Smina a la variable PATH")
    print("\n   O ejecutar:")
    print("   set PATH=%PATH%;C:\\ruta\\a\\smina")
    
    return False

def check_openbabel():
    """Verifica OpenBabel"""
    print("\n3. Verificando OpenBabel...")
    
    try:
        result = subprocess.run(
            ["obabel", "--version"],
            capture_output=True,
            text=True,
            timeout=5
        )
        if result.returncode == 0:
            version = result.stdout.split()[2] if len(result.stdout.split()) > 2 else "?"
            print(f"   ✓ OpenBabel {version} encontrado")
            return True
    except FileNotFoundError:
        pass
    except Exception:
        pass
    
    print("   ⚠ OpenBabel no encontrado (opcional)")
    print("      La preparación de ligandos se omitirá")
    return False

def check_python_packages():
    """Verifica paquetes Python opcionales"""
    print("\n4. Verificando paquetes Python opcionales...")
    
    packages = {
        'pandas': False,
        'matplotlib': False,
        'seaborn': False,
        'numpy': False
    }
    
    for package in packages:
        try:
            __import__(package)
            packages[package] = True
            print(f"   ✓ {package}")
        except ImportError:
            print(f"   ⚠ {package} no instalado (opcional)")
    
    if not any(packages.values()):
        print("\n   💡 Para instalar (opcional):")
        print("   pip install pandas matplotlib seaborn numpy")
    
    return packages

def check_files():
    """Verifica archivos de entrada"""
    print("\n5. Verificando archivos de entrada...")
    
    base_dir = Path.cwd()
    required_files = {
        '9FMM_protein.pdb': 'Receptor',
        'ligand_A1IDX_fragment2.pdb': 'Ligando Fragment2'
    }
    
    all_found = True
    for filename, description in required_files.items():
        filepath = base_dir / filename
        if filepath.exists():
            size = filepath.stat().st_size
            print(f"   ✓ {description}: {filename} ({size} bytes)")
        else:
            print(f"   ✗ {description}: {filename} NO ENCONTRADO")
            all_found = False
    
    if not all_found:
        print(f"\n   ⚠ Archivos faltantes en: {base_dir}")
    
    return all_found

def check_disk_space():
    """Verifica espacio en disco"""
    print("\n6. Verificando espacio en disco...")
    
    try:
        import shutil
        stat = shutil.disk_usage(Path.cwd())
        free_gb = stat.free / (1024**3)
        
        if free_gb >= 10:
            print(f"   ✓ Espacio libre: {free_gb:.1f} GB")
            return True
        else:
            print(f"   ⚠ Espacio libre: {free_gb:.1f} GB (recomendado: 10+ GB)")
            return False
    except Exception as e:
        print(f"   ⚠ No se pudo verificar espacio: {e}")
        return True

def main():
    """Main verification"""
    print("""
╔═══════════════════════════════════════════════════════════════════════════════╗
║                        VERIFICACIÓN DE ENTORNO                                ║
╚═══════════════════════════════════════════════════════════════════════════════╝
    """)
    
    print(f"Directorio actual: {Path.cwd()}")
    
    results = {
        'python': check_python(),
        'smina': check_smina(),
        'openbabel': check_openbabel(),
        'packages': check_python_packages(),
        'files': check_files(),
        'disk': check_disk_space()
    }
    
    print("\n" + "="*70)
    print("RESUMEN")
    print("="*70)
    
    # Critical requirements
    critical_ok = results['python'] and results['smina'] and results['files']
    
    if critical_ok:
        print("\n✅ REQUERIMIENTOS CRÍTICOS CUMPLIDOS")
        print("\nPuedes ejecutar:")
        print("  python master_pipeline_fixed.py")
    else:
        print("\n❌ FALTAN REQUERIMIENTOS CRÍTICOS:")
        if not results['python']:
            print("  - Actualizar Python a 3.7+")
        if not results['smina']:
            print("  - Instalar y configurar Smina")
        if not results['files']:
            print("  - Verificar archivos de entrada en directorio actual")
    
    # Optional features
    optional_features = []
    if not results['openbabel']:
        optional_features.append("Preparación de ligandos (OpenBabel)")
    
    any_package = any(results['packages'].values())
    if not any_package:
        optional_features.append("Análisis avanzado (pandas, matplotlib)")
    
    if optional_features:
        print("\n⚠️  CARACTERÍSTICAS OPCIONALES DESHABILITADAS:")
        for feature in optional_features:
            print(f"  - {feature}")
        print("\nEl pipeline funcionará pero sin estas características.")
    
    print("\n" + "="*70)
    
    if critical_ok:
        print("\n🚀 TODO LISTO PARA EMPEZAR!")
        print("\nEjecuta: python master_pipeline_fixed.py")
    else:
        print("\n⚠️  CONFIGURA LOS REQUERIMIENTOS CRÍTICOS ANTES DE CONTINUAR")

if __name__ == "__main__":
    main()