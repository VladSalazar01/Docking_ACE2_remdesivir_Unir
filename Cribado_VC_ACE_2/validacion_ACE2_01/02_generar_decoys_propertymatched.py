#!/usr/bin/env python3
"""
02_generar_decoys_propertymatched.py
====================================
Genera decoys property-matched para validación retrospectiva de docking.

METODOLOGÍA DUD-E ADAPTADA:
---------------------------
Los decoys son compuestos que:
1. Tienen propiedades fisicoquímicas SIMILARES a los activos
2. Tienen topología 2D DIFERENTE (no análogos estructurales)
3. No tienen actividad conocida contra la diana

Esto evita sesgos de propiedades que inflarían artificialmente las métricas.

CRITERIOS DE PROPERTY-MATCHING:
- Peso molecular: ±25%
- LogP: ±1.5
- Donadores H: ±2
- Aceptores H: ±3
- Carga neta: igual

CRITERIO DE DISIMILITUD TOPOLÓGICA:
- Coeficiente de Tanimoto < 0.35 (fingerprints Morgan)
- Verificado contra TODOS los activos del dataset

FUENTE DE SCAFFOLDS:
- Biblioteca curada de estructuras drug-like
- Sin actividad reportada contra ACE2 o proteínas relacionadas
- Diversidad estructural (antihistamínicos, heterociclos, ácidos, etc.)

REFERENCIAS:
- Mysinger MM et al. (2012) J Med Chem 55:6582-6594 - DUD-E methodology
- Empereur-Mot C et al. (2015) J Chem Inf Model 55:2066-2075 - Decoy generation

Autor: Script generado para TFM - Validación Gnina ACE2
Fecha: 2025-01-25
"""

import json
import os
import random
from pathlib import Path
from typing import Dict, List, Tuple

# Intentar importar dependencias
try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Crippen
except ImportError:
    print("Instalando RDKit...")
    os.system("pip install rdkit --break-system-packages -q")
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Crippen


# ============================================================================
# CONFIGURACIÓN
# ============================================================================
INPUT_DIR = Path("dataset_ace2/validacion")
OUTPUT_DIR = Path("dataset_ace2/validacion")

# Parámetros de matching (basados en DUD-E)
TOLERANCIA_MW = 0.25      # ±25% peso molecular
TOLERANCIA_LOGP = 1.5     # ±1.5 unidades LogP
TOLERANCIA_HBD = 2        # ±2 donadores H
TOLERANCIA_HBA = 3        # ±3 aceptores H
TANIMOTO_THRESHOLD = 0.35 # Máxima similaridad permitida

# Decoys objetivo por activo
DECOYS_POR_ACTIVO = 50


# ============================================================================
# BIBLIOTECA DE SCAFFOLDS DECOY
# ============================================================================
# Compuestos drug-like sin actividad conocida contra ACE2
# Organizados por familia química para diversidad

SCAFFOLDS_DECOY = {
    "antihistaminicos": [
        "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12",      # Cloroquina-like
        "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21",       # Fenotiazina
        "CC(C)Cc1ccc(C(C)C(=O)O)cc1",             # Ibuprofeno
        "Cc1ccccc1NC(=O)C(C)O",                   # Paracetamol-like
    ],
    
    "benzodiazepinas": [
        "c1ccc2c(c1)[nH]c1ccccc1c2=O",            # Quinazolinona
        "O=C1NC(=O)c2ccccc2N1",                   # Quinazolin-2,4-diona
        "c1ccc2[nH]ccc2c1",                       # Indol
        "c1ccc2ncccc2c1",                         # Quinolina
    ],
    
    "antiinflamatorios": [
        "Cc1ccc2c(C)c(C(=O)O)c(C)nc2c1",          # Derivado quinolina
        "O=c1[nH]c2ccccc2[nH]1",                  # Benzimidazolinona
        "O=C1CCc2ccccc2O1",                       # Croman-2-ona
        "c1ccc2occc2c1",                          # Benzofurano
    ],
    
    "heterociclos": [
        "c1cnc2[nH]ccc2n1",                       # Imidazo[1,2-a]pirimidina
        "c1cc2ccccc2[nH]1",                       # Indol alternativo
        "O=C1NCc2ccccc2N1",                       # Quinazolin-4(3H)-ona
        "c1ccc2sccc2c1",                          # Benzotiofeno
        "c1cnc2ccccc2n1",                         # Quinazolina
    ],
    
    "sulfonamidas": [
        "Nc1ccc(S(=O)(=O)Nc2ncccn2)cc1",          # Sulfadiazina-like
        "Cc1ncc(C(=O)O)cn1",                      # Piridina-carboxílico
    ],
    
    "nucleosidos": [
        "Nc1ncnc2[nH]cnc12",                      # Adenina
        "Nc1ccn(C)c(=O)n1",                       # Citosina N-metil
        "O=c1[nH]c(=O)c2[nH]cnc2[nH]1",           # Xantina
        "Cc1c[nH]c(=O)[nH]c1=O",                  # Timina
    ],
    
    "terpenos": [
        "CC1=CCC(C(C)(C)O)CC1",                   # Linalool-like
        "CC(=O)CC/C=C(/C)CCC=C(C)C",              # Geranilacetona-like
    ],
    
    "acidos_organicos": [
        "CCCCCCCCCCCC(=O)O",                      # Ácido láurico
        "CCCCC(=O)O",                             # Ácido valérico
        "CC(O)C(=O)O",                            # Ácido láctico
        "OC(=O)CC(O)(CC(=O)O)C(=O)O",             # Ácido cítrico
    ],
}

# Sustituciones para generar variantes
SUSTITUCIONES = [
    ("Cl", "F"), ("Cl", "Br"), ("F", "Cl"), ("F", "H"),
    ("O", "S"), ("N", "O"),
    ("C(=O)O", "C(=O)N"), ("C(=O)O", "CN"),
    ("OC", "SC"), ("NC", "OC"),
    ("CC", "CCC"), ("CCC", "CC"),
]


def calcular_propiedades(mol) -> Dict:
    """
    Calcula propiedades fisicoquímicas para property-matching.
    
    Propiedades según metodología DUD-E:
    - Peso molecular (MW)
    - LogP calculado (cLogP)
    - Donadores de enlace de hidrógeno (HBD)
    - Aceptores de enlace de hidrógeno (HBA)
    - Carga formal neta
    - Enlaces rotables
    """
    if mol is None:
        return None
    
    try:
        props = {
            "mw": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "rotatable": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        }
        return props
    except Exception:
        return None


def propiedades_similares(props_activo: Dict, props_decoy: Dict) -> bool:
    """
    Verifica si las propiedades son similares dentro de las tolerancias DUD-E.
    
    Criterios:
    - MW: ±25%
    - LogP: ±1.5
    - HBD: ±2
    - HBA: ±3
    - Carga: debe ser igual
    """
    if props_activo is None or props_decoy is None:
        return False
    
    # Peso molecular: ±25%
    mw_diff = abs(props_decoy["mw"] - props_activo["mw"])
    if mw_diff > props_activo["mw"] * TOLERANCIA_MW:
        return False
    
    # LogP: ±1.5
    if abs(props_decoy["logp"] - props_activo["logp"]) > TOLERANCIA_LOGP:
        return False
    
    # HBD: ±2
    if abs(props_decoy["hbd"] - props_activo["hbd"]) > TOLERANCIA_HBD:
        return False
    
    # HBA: ±3
    if abs(props_decoy["hba"] - props_activo["hba"]) > TOLERANCIA_HBA:
        return False
    
    # Carga: debe ser igual
    if props_decoy["charge"] != props_activo["charge"]:
        return False
    
    return True


def calcular_fingerprint(mol):
    """
    Calcula fingerprint Morgan (ECFP4-like) para comparación topológica.
    
    Parámetros:
    - Radio: 2 (equivalente a ECFP4)
    - nBits: 2048
    """
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)


def topologia_diferente(fps_activos: List, fp_decoy) -> bool:
    """
    Verifica que el decoy sea topológicamente diferente a TODOS los activos.
    
    Criterio: Tanimoto < 0.35 contra todos los activos
    
    Esto asegura que los decoys no sean análogos estructurales que 
    podrían tener actividad similar.
    """
    for fp_activo in fps_activos:
        similaridad = DataStructs.TanimotoSimilarity(fp_activo, fp_decoy)
        if similaridad >= TANIMOTO_THRESHOLD:
            return False
    
    return True


def generar_variantes(smiles_base: str, max_variantes: int = 20) -> List[str]:
    """
    Genera variantes estructurales mediante sustituciones bioisostéricas.
    
    Esto expande la diversidad del pool de decoys candidatos.
    """
    variantes = []
    mol = Chem.MolFromSmiles(smiles_base)
    
    if mol is None:
        return variantes
    
    for original, reemplazo in SUSTITUCIONES:
        if original in smiles_base:
            nuevo_smiles = smiles_base.replace(original, reemplazo, 1)
            mol_nuevo = Chem.MolFromSmiles(nuevo_smiles)
            
            if mol_nuevo is not None:
                try:
                    smiles_can = Chem.MolToSmiles(mol_nuevo)
                    if smiles_can not in variantes:
                        variantes.append(smiles_can)
                except Exception:
                    pass
    
    return variantes[:max_variantes]


def construir_pool_decoys() -> List[Dict]:
    """
    Construye el pool de candidatos a decoys desde scaffolds y variantes.
    """
    print("\n" + "-" * 50)
    print("CONSTRUCCIÓN DE POOL DE DECOYS CANDIDATOS")
    print("-" * 50)
    
    pool = []
    
    for familia, scaffolds in SCAFFOLDS_DECOY.items():
        count_familia = 0
        
        for scaffold in scaffolds:
            mol = Chem.MolFromSmiles(scaffold)
            if mol is None:
                continue
            
            props = calcular_propiedades(mol)
            if props and props["mw"] < 600:  # Filtrar MW razonable
                pool.append({
                    "smiles": Chem.MolToSmiles(mol),
                    "mol": mol,
                    "props": props,
                    "familia": familia
                })
                count_familia += 1
            
            # Generar variantes
            for variante in generar_variantes(scaffold, 20):
                mol_var = Chem.MolFromSmiles(variante)
                if mol_var is not None:
                    props_var = calcular_propiedades(mol_var)
                    if props_var and props_var["mw"] < 600:
                        pool.append({
                            "smiles": variante,
                            "mol": mol_var,
                            "props": props_var,
                            "familia": familia
                        })
                        count_familia += 1
        
        print(f"  {familia}: {count_familia} candidatos")
    
    print(f"\nTotal pool: {len(pool)} candidatos")
    return pool


def cargar_activos() -> List[Dict]:
    """
    Carga activos desde archivo SMI y calcula propiedades/fingerprints.
    """
    activos_file = INPUT_DIR / "activos.smi"
    
    if not activos_file.exists():
        print(f"ERROR: No se encuentra {activos_file}")
        print("Ejecuta primero: python 01_curar_activos_literatura.py")
        return []
    
    activos = []
    
    with open(activos_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                smiles = parts[0]
                nombre = parts[1]
                
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    activos.append({
                        "smiles": smiles,
                        "nombre": nombre,
                        "mol": mol,
                        "props": calcular_propiedades(mol),
                        "fp": calcular_fingerprint(mol)
                    })
    
    return activos


def seleccionar_decoys(activos: List[Dict], pool: List[Dict]) -> List[Dict]:
    """
    Selecciona decoys property-matched y topológicamente diferentes.
    
    Algoritmo:
    1. Para cada activo, iterar sobre pool aleatorizado
    2. Verificar property-matching
    3. Verificar disimilitud topológica vs TODOS los activos
    4. Evitar duplicados
    5. Continuar hasta alcanzar cuota o agotar pool
    """
    print("\n" + "-" * 50)
    print("SELECCIÓN DE DECOYS PROPERTY-MATCHED")
    print("-" * 50)
    
    # Fingerprints de todos los activos
    fps_activos = [a["fp"] for a in activos]
    
    decoys_seleccionados = []
    smiles_usados = set()
    
    # Aleatorizar pool
    random.shuffle(pool)
    
    for activo in activos:
        count = 0
        
        for candidato in pool:
            if count >= DECOYS_POR_ACTIVO:
                break
            
            # Evitar duplicados
            if candidato["smiles"] in smiles_usados:
                continue
            
            # Verificar property-matching
            if not propiedades_similares(activo["props"], candidato["props"]):
                continue
            
            # Verificar topología diferente
            fp_candidato = calcular_fingerprint(candidato["mol"])
            if not topologia_diferente(fps_activos, fp_candidato):
                continue
            
            # Seleccionar
            decoys_seleccionados.append({
                "smiles": candidato["smiles"],
                "familia": candidato["familia"],
                "matched_to": activo["nombre"]
            })
            smiles_usados.add(candidato["smiles"])
            count += 1
        
        print(f"  {activo['nombre']}: {count} decoys seleccionados")
    
    return decoys_seleccionados


def guardar_decoys(decoys: List[Dict], activos: List[Dict]):
    """
    Guarda decoys en formato SMI y genera resumen estadístico.
    """
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Archivo SMI
    decoys_file = OUTPUT_DIR / "decoys.smi"
    with open(decoys_file, "w") as f:
        for i, decoy in enumerate(decoys):
            f.write(f"{decoy['smiles']} decoy_{i:04d}\n")
    
    print(f"\nArchivo guardado: {decoys_file}")
    
    # Resumen estadístico
    resumen_file = OUTPUT_DIR / "resumen_dataset.txt"
    with open(resumen_file, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write("DATASET DE VALIDACIÓN ACE2\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Activos: {len(activos)}\n")
        f.write(f"Decoys:  {len(decoys)}\n")
        f.write(f"Ratio:   1:{len(decoys) // max(len(activos), 1)}\n\n")
        
        f.write("-" * 40 + "\n")
        f.write("CRITERIOS DE SELECCIÓN DE DECOYS:\n")
        f.write("-" * 40 + "\n")
        f.write(f"MW tolerance:     ±{TOLERANCIA_MW*100:.0f}%\n")
        f.write(f"LogP tolerance:   ±{TOLERANCIA_LOGP}\n")
        f.write(f"HBD tolerance:    ±{TOLERANCIA_HBD}\n")
        f.write(f"HBA tolerance:    ±{TOLERANCIA_HBA}\n")
        f.write(f"Tanimoto cutoff:  <{TANIMOTO_THRESHOLD}\n\n")
        
        f.write("-" * 40 + "\n")
        f.write("DISTRIBUCIÓN POR FAMILIA:\n")
        f.write("-" * 40 + "\n")
        familias = {}
        for d in decoys:
            fam = d["familia"]
            familias[fam] = familias.get(fam, 0) + 1
        for fam, count in sorted(familias.items(), key=lambda x: -x[1]):
            f.write(f"  {fam}: {count}\n")
    
    print(f"Resumen guardado: {resumen_file}")


def main():
    """Función principal."""
    print("=" * 70)
    print("GENERACIÓN DE DECOYS PROPERTY-MATCHED")
    print("=" * 70)
    
    # 1. Cargar activos
    print("\nCargando activos...")
    activos = cargar_activos()
    
    if not activos:
        print("ERROR: No se pudieron cargar activos")
        return
    
    print(f"Activos cargados: {len(activos)}")
    
    # 2. Construir pool de decoys
    pool = construir_pool_decoys()
    
    # 3. Seleccionar decoys
    decoys = seleccionar_decoys(activos, pool)
    
    # 4. Guardar resultados
    guardar_decoys(decoys, activos)
    
    # 5. Resumen final
    print("\n" + "=" * 70)
    print("RESUMEN DEL DATASET DE VALIDACIÓN ACE2")
    print("=" * 70)
    print(f"Activos: {len(activos)}")
    print(f"Decoys:  {len(decoys)}")
    print(f"Ratio:   1:{len(decoys) // max(len(activos), 1)}")
    
    if len(decoys) < len(activos) * 10:
        print("\n⚠️  NOTA: Ratio bajo (<1:10)")
        print("   Para mejor validación, considerar:")
        print("   - DUD-E Server (https://dude.docking.org/generate)")
        print("   - Descargar ZINC leads para mayor diversidad")
    
    print("\n" + "=" * 70)
    print("SIGUIENTE PASO: EJECUTAR VALIDACIÓN")
    print("=" * 70)
    print("""
Ejecutar docking de validación:

  python run_gnina_validacion_ace2.py

O calcular métricas si ya se ejecutó docking:

  python calcular_metricas.py
""")


if __name__ == "__main__":
    main()
