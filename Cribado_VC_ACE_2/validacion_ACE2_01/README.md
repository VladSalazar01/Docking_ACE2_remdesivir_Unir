# Validación Retrospectiva Gnina - ACE2

## Resumen del Proyecto

Este directorio contiene los scripts y dataset para la validación retrospectiva del protocolo Gnina en la diana ACE2 (Angiotensin-converting enzyme 2), como parte del TFM sobre cribado virtual de moduladores ACE2.

### Contexto

Dado que ACE2 tiene **muy pocos inhibidores small-molecule documentados**, esta validación complementa la validación cruzada ya realizada con ACE1 (homólogo con ~37% identidad de secuencia).

**Inhibidores ACE2 documentados:**
- **MLN-4760**: IC50 = 0.44 nM (único inhibidor potente robusto)
- **Flavonoides**: Quercetina (IC50 ~4.5 μM), Rutina (~12 μM), etc.
- **Péptidos**: DX600, péptidos fosfínicos (no aptos para docking small-molecule)

## Estructura de Archivos

```
validacion_ace2/
├── README.md                           # Este archivo
├── obtener_activos_chembl.py          # Script para buscar activos en ChEMBL
├── generar_decoys_local.py            # Generador de decoys property-matched
├── run_gnina_validacion_ace2.py       # Script principal de validación
├── calcular_metricas.py               # Cálculo de métricas post-docking
└── dataset_ace2/
    ├── activos_ace2.smi               # Activos básicos (4 compuestos)
    ├── activos_ace2_expandido.smi     # Dataset expandido con literatura
    ├── activos_ace2_metadata.json     # Metadatos JSON
    ├── resumen_activos.txt            # Resumen estadístico
    └── validacion/
        ├── activos.smi                # Activos finales procesados (19)
        └── decoys.smi                 # Decoys property-matched (41)
```

## Dataset de Validación

### Composición
- **Activos**: 19 compuestos
  - 2 muy activos (IC50 ≤100 nM): MLN-4760, F-MLN-4760
  - 12 moderados (1-100 μM): Flavonoides y ácidos fenólicos
  - 5 análogos estructurales de MLN-4760
  
- **Decoys**: 41 compuestos property-matched
  - Scaffolds drug-like sin actividad ACE2 conocida
  - Matching: MW (±25%), LogP (±1.5), HBD (±2), HBA (±3)
  - Tanimoto < 0.35 vs todos los activos

- **Ratio**: 1:2 (limitado por disponibilidad de scaffolds)

### Nota sobre el Dataset
Para una validación más robusta, se recomienda:
1. Usar DUD-E Server (https://dude.docking.org/generate) para generar 50 decoys/activo
2. O descargar ZINC leads para mayor diversidad

## Instrucciones de Ejecución

### 1. Preparar Receptor ACE2

Usar el receptor PDB 9FMM ya preparado del benchmark previo:
```bash
# Si no está preparado:
cd G:\Cribado_VC_ACE_2\validacion_ACE2_01
# Copiar 9FMM_receptor.pdbqt a la carpeta de validación
```

### 2. Copiar Dataset a Windows

```batch
xcopy /E /I validacion_ace2 G:\Cribado_VC_ACE_2\validacion_ACE2_01
```

### 3. Ejecutar Docking

**Opción A - Script batch (recomendado):**
```batch
cd G:\Cribado_VC_ACE_2\validacion_ACE2_01
.\ejecutar_docking.bat
```

**Opción B - Comando individual:**
```batch
docker run --rm --gpus all -v "%cd%":/data gnina/gnina gnina ^
  --receptor /data/9FMM_receptor.pdbqt ^
  --ligand /data/ligandos_pdbqt/MLN4760.pdbqt ^
  --out /data/salidas/MLN4760_out.pdbqt ^
  --center_x 42.0 --center_y 7.0 --center_z 23.0 ^
  --size_x 26 --size_y 26 --size_z 26 ^
  --exhaustiveness 64 ^
  --cnn_scoring rescore ^
  --num_modes 1
```

### 4. Calcular Métricas

```bash
python calcular_metricas.py
```

## Parámetros de Docking

| Parámetro | Valor | Notas |
|-----------|-------|-------|
| Receptor | PDB 9FMM | ACE2 humana + F-MLN-4760 |
| Centro X | 42.0 Å | Desde ligando cristalográfico |
| Centro Y | 7.0 Å | |
| Centro Z | 23.0 Å | |
| Tamaño | 26×26×26 Å | Cubre sitio activo completo |
| Exhaustiveness | 64 | Optimizado en benchmark previo |
| CNN Scoring | rescore | Mejor para screening |
| Num Modes | 1 | Solo mejor pose para validación |

## Umbrales de Aceptación

| Métrica | Umbral Mínimo | Descripción |
|---------|---------------|-------------|
| AUC-ROC | ≥ 0.70 | Capacidad discriminatoria global |
| EF 1% | ≥ 10 | Enriquecimiento en top 1% |
| EF 5% | ≥ 5 | Enriquecimiento en top 5% |
| EF 10% | ≥ 3 | Enriquecimiento en top 10% |
| BEDROC | ≥ 0.30 | Enriquecimiento temprano (α=20) |

**Criterio de validación**: Mínimo 3 de 5 métricas cumplidas

## Comparación con Validación Cruzada ACE1

| Métrica | ACE1 (Cross-val) | ACE2 (Esperado) |
|---------|------------------|-----------------|
| AUC-ROC | 0.641 | 0.65-0.75 |
| EF 1% | 6.82 | 5-10 |
| EF 5% | 3.64 | 3-6 |
| BEDROC | 0.191 | 0.20-0.35 |

## Justificación Metodológica

1. **Dataset pequeño**: ACE2 tiene pocos inhibidores small-molecule publicados, lo cual es precisamente por qué la validación cruzada con ACE1 (diana relacionada con dataset DUD-E disponible) es importante.

2. **Inhibidor de referencia**: MLN-4760 (IC50 = 0.44 nM) es el gold standard para validación. Su presencia en el top de ranking es crucial.

3. **Complementariedad con ACE1**: La homología estructural (~37%) justifica usar ambas validaciones como evidencia de la robustez del protocolo.

## Tiempo Estimado

- Preparación de ligandos: ~5 minutos
- Docking (60 compuestos, 2×RTX3060): ~30 minutos
- Cálculo de métricas: ~1 minuto

**Total**: ~40 minutos

## Archivos de Salida

```
resultados_validacion_ace2/
├── ligandos_pdbqt/        # Ligandos preparados
├── salidas/               # Poses de docking
├── resultados/            # Logs individuales
└── metricas_finales/
    ├── reporte_validacion_ace2.txt
    ├── metricas_validacion_ace2.json
    ├── curva_roc_ace2.png
    └── metricas_enriquecimiento_ace2.png
```

## Referencias

1. Dales NA et al. (2002) JACS 124:11852 - Diseño de MLN-4760
2. Towler P et al. (2004) JBC 279:17996 - Estructura cristalográfica ACE2
3. Mysinger MM et al. (2012) J Med Chem 55:6582 - Metodología DUD-E
4. Ragoza M et al. (2017) J Chem Inf Model 57:942 - Gnina

---
*Generado para TFM - Validación Gnina ACE2*
*Fecha: 2025-01-25*
