#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
diagnostico_atomos.py
Diagnóstico del desajuste atómico entre referencia PDB y output PDBQT de Vina.
Identifica si el problema de RMSD alto (~5 Å) es por reordenamiento de átomos.
"""

from pathlib import Path
import numpy as np

CARPETA = Path(r"G:\Autodock_vina_benchmark")

def leer_atomos(ruta, formato="pdb"):
    """Lee átomos de PDB o PDBQT. Retorna lista de dicts."""
    atomos = []
    with open(ruta, "r") as f:
        en_primer_modelo = True
        for linea in f:
            if linea.startswith("ENDMDL"):
                break  # solo primer modelo en PDBQT multi-pose
            if linea.startswith("MODEL") and len(atomos) > 0:
                break
            if linea.startswith("HETATM") or linea.startswith("ATOM"):
                nombre = linea[12:16].strip()
                resname = linea[17:20].strip()
                x = float(linea[30:38])
                y = float(linea[38:46])
                z = float(linea[46:54])
                # Elemento
                if formato == "pdbqt" and len(linea) > 77:
                    elemento = linea[77:79].strip()
                elif len(linea) > 76:
                    elemento = linea[76:78].strip()
                else:
                    elemento = nombre[0]
                atomos.append({
                    "nombre": nombre,
                    "resname": resname,
                    "x": x, "y": y, "z": z,
                    "elemento": elemento,
                    "linea_raw": linea.rstrip()
                })
    return atomos

def es_hidrogeno(atomo):
    """Determina si un átomo es hidrógeno."""
    el = atomo["elemento"].upper()
    nom = atomo["nombre"].upper()
    return el in ("H", "HD") or (len(el) == 0 and nom.startswith("H"))

def filtrar_pesados(atomos):
    """Retorna solo átomos pesados (no H)."""
    return [a for a in atomos if not es_hidrogeno(a)]

def centro_masa(atomos):
    """Calcula centro de masa."""
    coords = np.array([[a["x"], a["y"], a["z"]] for a in atomos])
    return coords.mean(axis=0)

def kabsch_rmsd(P, Q):
    """RMSD con alineamiento Kabsch."""
    n = min(len(P), len(Q))
    P, Q = P[:n], Q[:n]
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign = np.diag([1, 1, np.sign(d)])
    R = Vt.T @ sign @ U.T
    Prot = (R @ Pc.T).T
    return np.sqrt(((Prot - Qc)**2).sum() / n)

def rmsd_directo(P, Q):
    """RMSD sin alineamiento (solo para ver magnitud bruta)."""
    n = min(len(P), len(Q))
    return np.sqrt(((P[:n] - Q[:n])**2).sum() / n)

def intentar_match_por_nombre(ref_pesados, dock_pesados):
    """Intenta emparejar átomos por nombre y calcula RMSD."""
    ref_dict = {}
    for a in ref_pesados:
        key = a["nombre"]
        if key not in ref_dict:
            ref_dict[key] = a
    
    pareados_ref = []
    pareados_dock = []
    no_encontrados = []
    
    for a in dock_pesados:
        key = a["nombre"]
        if key in ref_dict:
            r = ref_dict[key]
            pareados_ref.append([r["x"], r["y"], r["z"]])
            pareados_dock.append([a["x"], a["y"], a["z"]])
        else:
            no_encontrados.append(key)
    
    return np.array(pareados_ref), np.array(pareados_dock), no_encontrados

def diagnosticar_sitio(nombre_sitio, ref_pdb, lig_pdbqt, dock_pdbqt):
    """Diagnóstico completo para un sitio."""
    print(f"\n{'='*70}")
    print(f"  DIAGNÓSTICO: {nombre_sitio}")
    print(f"{'='*70}")
    
    # --- Leer archivos ---
    ref_atomos = leer_atomos(ref_pdb, "pdb") if ref_pdb.exists() else []
    lig_atomos = leer_atomos(lig_pdbqt, "pdbqt") if lig_pdbqt.exists() else []
    dock_atomos = leer_atomos(dock_pdbqt, "pdbqt") if dock_pdbqt.exists() else []
    
    if not ref_atomos:
        print(f"  [ERROR] Referencia no encontrada: {ref_pdb}")
        return
    if not dock_atomos:
        print(f"  [ERROR] Pose docked no encontrada: {dock_pdbqt}")
        return
    
    ref_pesados = filtrar_pesados(ref_atomos)
    lig_pesados = filtrar_pesados(lig_atomos) if lig_atomos else []
    dock_pesados = filtrar_pesados(dock_atomos)
    
    # --- Conteos ---
    print(f"\n  1. CONTEO DE ÁTOMOS")
    print(f"  {'Archivo':<40} {'Total':<8} {'Pesados':<8} {'H':<8}")
    print(f"  {'-'*64}")
    print(f"  {'Referencia PDB':<40} {len(ref_atomos):<8} {len(ref_pesados):<8} {len(ref_atomos)-len(ref_pesados):<8}")
    if lig_atomos:
        print(f"  {'Ligando PDBQT (input)':<40} {len(lig_atomos):<8} {len(lig_pesados):<8} {len(lig_atomos)-len(lig_pesados):<8}")
    print(f"  {'Pose docked PDBQT':<40} {len(dock_atomos):<8} {len(dock_pesados):<8} {len(dock_atomos)-len(dock_pesados):<8}")
    
    if len(ref_pesados) != len(dock_pesados):
        print(f"\n  [AVISO]  DESAJUSTE: ref={len(ref_pesados)} vs dock={len(dock_pesados)} átomos pesados")
    else:
        print(f"\n  [OK] Mismo número de átomos pesados: {len(ref_pesados)}")
    
    # --- Nombres de átomos ---
    print(f"\n  2. ORDEN DE ÁTOMOS (pesados)")
    print(f"  {'#':<4} {'Ref (PDB)':<12} {'Dock (PDBQT)':<12} {'¿Coincide?'}")
    print(f"  {'-'*40}")
    n_show = max(len(ref_pesados), len(dock_pesados))
    coinciden = 0
    for i in range(min(n_show, 35)):
        r_name = ref_pesados[i]["nombre"] if i < len(ref_pesados) else "---"
        d_name = dock_pesados[i]["nombre"] if i < len(dock_pesados) else "---"
        match = "[OK]" if r_name == d_name else "[X]"
        if r_name == d_name:
            coinciden += 1
        print(f"  {i+1:<4} {r_name:<12} {d_name:<12} {match}")
    if n_show > 35:
        print(f"  ... ({n_show - 35} átomos más)")
    
    total_comparables = min(len(ref_pesados), len(dock_pesados))
    print(f"\n  Coincidencia de nombres: {coinciden}/{total_comparables} "
          f"({100*coinciden/total_comparables:.0f}%)" if total_comparables > 0 else "")
    
    # --- Centros de masa ---
    print(f"\n  3. CENTROS DE MASA")
    cm_ref = centro_masa(ref_pesados)
    cm_dock = centro_masa(dock_pesados)
    dist_cm = np.linalg.norm(cm_ref - cm_dock)
    print(f"  Referencia:  ({cm_ref[0]:.3f}, {cm_ref[1]:.3f}, {cm_ref[2]:.3f})")
    print(f"  Docked:      ({cm_dock[0]:.3f}, {cm_dock[1]:.3f}, {cm_dock[2]:.3f})")
    print(f"  Distancia:   {dist_cm:.3f} Å")
    
    if dist_cm > 10:
        print(f"  [AVISO]  ALERTA: Centro de masa >10 Å de diferencia — ¿sitio equivocado?")
    elif dist_cm > 5:
        print(f"  [AVISO]  AVISO: Centro de masa >5 Å — posible problema de posicionamiento")
    
    # --- RMSD por diferentes métodos ---
    print(f"\n  4. RMSD POR DIFERENTES MÉTODOS")
    
    coords_ref = np.array([[a["x"], a["y"], a["z"]] for a in ref_pesados])
    coords_dock = np.array([[a["x"], a["y"], a["z"]] for a in dock_pesados])
    
    # 4a. RMSD directo (sin alinear, orden original)
    if len(coords_ref) == len(coords_dock):
        rmsd_dir = rmsd_directo(coords_ref, coords_dock)
        print(f"  a) Directo (sin alinear):       {rmsd_dir:.3f} Å")
    
    # 4b. Kabsch (orden original)
    rmsd_k = kabsch_rmsd(coords_ref, coords_dock)
    print(f"  b) Kabsch (orden original):     {rmsd_k:.3f} Å  ← lo que reporta el script")
    
    # 4c. Match por nombre de átomo + Kabsch
    ref_matched, dock_matched, no_encontrados = intentar_match_por_nombre(ref_pesados, dock_pesados)
    if len(ref_matched) > 0:
        rmsd_name = kabsch_rmsd(ref_matched, dock_matched)
        print(f"  c) Kabsch (matched por nombre): {rmsd_name:.3f} Å  ({len(ref_matched)} pares)")
        if no_encontrados:
            print(f"     No emparejados: {no_encontrados[:10]}")
    
    # 4d. Sorted por elemento + coordenada (heurística)
    def sort_key(a):
        return (a["elemento"], round(a["x"], 1), round(a["y"], 1))
    
    ref_sorted = sorted(ref_pesados, key=sort_key)
    dock_sorted = sorted(dock_pesados, key=sort_key)
    n_sorted = min(len(ref_sorted), len(dock_sorted))
    coords_ref_s = np.array([[a["x"], a["y"], a["z"]] for a in ref_sorted[:n_sorted]])
    coords_dock_s = np.array([[a["x"], a["y"], a["z"]] for a in dock_sorted[:n_sorted]])
    rmsd_sorted = kabsch_rmsd(coords_ref_s, coords_dock_s)
    print(f"  d) Kabsch (sorted heurístico):  {rmsd_sorted:.3f} Å")
    
    # --- Diagnóstico ---
    print(f"\n  5. DIAGNÓSTICO")
    if coinciden < total_comparables * 0.5:
        print(f"  [!!] PROBLEMA CONFIRMADO: Orden de átomos DIFERENTE entre ref y docked")
        print(f"     obabel reordena átomos al convertir PDB→PDBQT")
        print(f"     Solución: usar match por nombre de átomo (método c)")
        if len(ref_matched) > 0 and rmsd_name < rmsd_k * 0.5:
            print(f"     → RMSD real estimado: {rmsd_name:.3f} Å (vs {rmsd_k:.3f} Å reportado)")
    elif dist_cm > 10:
        print(f"  [!!] PROBLEMA: Pose docked muy lejos de la referencia")
        print(f"     ¿La caja de docking contiene el sitio correcto?")
    elif rmsd_k > 3.0:
        print(f"  [--] RMSD alto pero átomos coinciden. Probable limitación de Vina para este sistema.")
    else:
        print(f"  [OK] Sin problemas evidentes detectados")


def main():
    print("=" * 70)
    print("  DIAGNÓSTICO DE ÁTOMOS - BENCHMARK VINA DUAL-SITE")
    print("=" * 70)
    
    # Listar archivos disponibles
    print(f"\n  Archivos en {CARPETA}:")
    for ext in ("*.pdb", "*.pdbqt"):
        for f in sorted(CARPETA.glob(ext)):
            size_kb = f.stat().st_size / 1024
            print(f"    {f.name:<45} {size_kb:.1f} KB")
    
    # Diagnosticar cada sitio
    for nombre_sitio in ["Fragment_1", "Fragment_2"]:
        ref_pdb = CARPETA / f"ligand_exp_{nombre_sitio}.pdb"
        lig_pdbqt = CARPETA / f"ligand_{nombre_sitio}.pdbqt"
        
        # Buscar primer output docked
        dock_pdbqt = CARPETA / f"vina_out_{nombre_sitio}_001.pdbqt"
        if not dock_pdbqt.exists():
            # Buscar en carpeta de resultados más reciente
            resultados = sorted(CARPETA.glob(f"resultados_dual_site_*"), reverse=True)
            for carpeta_res in resultados:
                candidato = carpeta_res / f"vina_out_{nombre_sitio}_001.pdbqt"
                if candidato.exists():
                    dock_pdbqt = candidato
                    break
        
        diagnosticar_sitio(nombre_sitio, ref_pdb, lig_pdbqt, dock_pdbqt)
    
    # Verificar si existe el ligando completo original
    print(f"\n{'='*70}")
    print(f"  ARCHIVOS ADICIONALES")
    print(f"{'='*70}")
    for nombre in ["ligand.pdbqt", "ligand_exp.pdb", "ligand_A1IDX.pdb", "ligand_A1IDX.pdbqt"]:
        ruta = CARPETA / nombre
        if ruta.exists():
            atomos = leer_atomos(ruta, "pdbqt" if nombre.endswith(".pdbqt") else "pdb")
            pesados = filtrar_pesados(atomos)
            cm = centro_masa(pesados) if pesados else np.zeros(3)
            print(f"  [OK] {nombre:<35} {len(atomos)} átomos ({len(pesados)} pesados), "
                  f"CM=({cm[0]:.1f}, {cm[1]:.1f}, {cm[2]:.1f})")
        else:
            print(f"  [X] {nombre:<35} NO EXISTE")


if __name__ == "__main__":
    main()