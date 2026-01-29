# Sección para TFM: Validación Retrospectiva ACE2

## Texto propuesto para sección de Metodología

### 3.X.X Validación retrospectiva específica para ACE2

La validación directa del protocolo Gnina sobre ACE2 se realizó mediante un dataset curado de inhibidores con actividad enzimática documentada. A diferencia de ACE1, que dispone de un benchmark DUD-E establecido, ACE2 carece de un dataset de validación estándar debido a la escasez de inhibidores small-molecule publicados. Por ello, se construyó un conjunto de validación *ad hoc* siguiendo los principios metodológicos de DUD-E.

**Dataset de activos**: Se compilaron 19 compuestos con actividad inhibitoria demostrada contra ACE2 humana, incluyendo MLN-4760 (IC50 = 0.44 nM) como referencia principal, su análogo difluorado F-MLN-4760 presente en la estructura cristalográfica 9FMM, y una serie de flavonoides (quercetina, rutina, luteolina, kaempferol) con actividades en el rango micromolar. Los criterios de inclusión fueron: IC50 o Ki ≤ 100 μM, peso molecular ≤ 600 Da, y disponibilidad de estructura química validada.

**Dataset de decoys**: Se generaron decoys property-matched siguiendo criterios adaptados de DUD-E: peso molecular (±25%), LogP calculado (±1.5), donores de hidrógeno (±2), aceptores de hidrógeno (±3), y similaridad topológica Tanimoto < 0.35 respecto a todos los activos. El dataset final comprendió 41 decoys seleccionados de una biblioteca de scaffolds farmacológicamente diversos sin actividad ACE2 conocida.

**Parámetros de docking**: Se emplearon los parámetros optimizados en el benchmark previo: centro de grid (42.0, 7.0, 23.0) Å derivado del ligando cristalográfico, tamaño de caja 26×26×26 Å, exhaustiveness 64, y CNN scoring en modo *rescore*.

---

## Tabla propuesta para sección de Resultados

**Tabla X. Métricas de validación retrospectiva del protocolo Gnina para ACE2**

| Métrica | Valor obtenido | Umbral aceptable | Evaluación |
|---------|----------------|------------------|------------|
| AUC-ROC | [PENDIENTE] | ≥ 0.70 | [CUMPLE/NO CUMPLE] |
| EF 1% | [PENDIENTE] | ≥ 10 | [CUMPLE/NO CUMPLE] |
| EF 5% | [PENDIENTE] | ≥ 5 | [CUMPLE/NO CUMPLE] |
| EF 10% | [PENDIENTE] | ≥ 3 | [CUMPLE/NO CUMPLE] |
| BEDROC (α=20) | [PENDIENTE] | ≥ 0.30 | [CUMPLE/NO CUMPLE] |

*Nota: Dataset compuesto por 19 activos y 41 decoys (ratio 1:2.2). Receptor: PDB 9FMM.*

---

## Tabla comparativa con validación cruzada ACE1

**Tabla Y. Comparación de métricas de validación entre dianas ACE1 y ACE2**

| Métrica | ACE1 (DUD-E) | ACE2 (Dataset ad hoc) | Diferencia |
|---------|--------------|----------------------|------------|
| n activos | 50 | 19 | -31 |
| n decoys | 1000 | 41 | -959 |
| Ratio | 1:20 | 1:2.2 | - |
| AUC-ROC | 0.641 | [PENDIENTE] | - |
| EF 1% | 6.82 | [PENDIENTE] | - |
| EF 5% | 3.64 | [PENDIENTE] | - |
| BEDROC | 0.191 | [PENDIENTE] | - |

---

## Texto propuesto para Discusión

### 4.X Validación del protocolo en ACE2

Los resultados de la validación retrospectiva específica para ACE2 [complementan/confirman] los hallazgos obtenidos en la validación cruzada con ACE1. La diferencia en el tamaño de los datasets refleja la realidad bibliográfica de estas dianas: mientras ACE1 cuenta con cientos de inhibidores documentados en ChEMBL, ACE2 presenta una escasez notable de moduladores small-molecule, siendo MLN-4760 prácticamente el único inhibidor potente (IC50 < 1 nM) disponible.

[Si los resultados son positivos:]
El posicionamiento de MLN-4760 en el top del ranking de scores CNN valida la capacidad del protocolo para identificar el inhibidor más potente conocido. Los valores de AUC-ROC y EF obtenidos, aunque [superiores/inferiores] a los de ACE1, demuestran que el protocolo discrimina adecuadamente entre activos y decoys en esta diana específica.

[Si los resultados son subóptimos:]
Las métricas de validación para ACE2, aunque por debajo de los umbrales óptimos establecidos, deben interpretarse considerando las limitaciones inherentes del dataset: tamaño reducido, alto peso de compuestos micromolares con menor especificidad de binding, y ratio decoys:activos subóptimo. La validación cruzada con ACE1 proporciona evidencia complementaria de la robustez general del protocolo Gnina.

### 4.X.1 Limitaciones de la validación ACE2

La principal limitación de esta validación radica en la escasez de inhibidores small-molecule documentados para ACE2. A diferencia de otras dianas farmacológicas establecidas, ACE2 emergió como diana terapéutica prominente principalmente tras la pandemia COVID-19, con mayor énfasis en moduladores de la interacción Spike-ACE2 que en inhibidores enzimáticos clásicos. Esto resulta en un dataset de validación limitado tanto en número como en diversidad estructural, sesgado hacia dos familias químicas: derivados de MLN-4760 (potentes) y flavonoides (moderados).

---

## Notas para completar después del docking

1. Ejecutar scripts de validación según README.md
2. Rellenar valores [PENDIENTE] con resultados reales
3. Cambiar [CUMPLE/NO CUMPLE] según umbrales
4. Ajustar texto de discusión según resultados obtenidos
5. Incluir figuras: curva ROC y gráfica de enriquecimiento
