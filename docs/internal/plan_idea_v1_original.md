# smrna\_nextflow — Especificación técnica v1.0

*Pipeline downstream de small RNA-seq sobre los resultados de nf-core/smrnaseq. Proyecto CRUK\_DDX41\_DHX34, batch\_2026\_04. Documento pensado para que lo implemente un junior o un agente sin ambigüedades.*

## 0. Cómo leer este documento

- **\[OBLIGATORIO\]**: debe cumplirse tal cual. **\[DEFAULT\]**: valor por defecto; solo se cambia vía parámetro. **\[VERIFICAR\]**: dato que no se ha podido comprobar desde aquí; hay que comprobarlo en el primer contacto con los datos y, si no coincide, parar y reportar (no improvisar). **\[OPCIONAL\]**: módulo o salida que puede quedar desactivado.
- Regla general: si algo no está especificado, **parar y preguntar**. No inventar columnas, umbrales ni rutas.
- Código, parámetros y nombres de fichero en inglés. Los datos de entrada son de solo lectura: ningún módulo escribe fuera de `outdir`.

## 1. Objetivo y alcance

**Contrato:** entrada = carpeta `results/` de nf-core/smrnaseq 2.4.x + `metadata.csv` + `contrasts.csv`; salida = carpeta `results_smrna/` con tablas, figuras e informe.

| Dentro del alcance | Fuera del alcance (de momento) |
| --- | --- |
| QC por muestra y flags | Repetir recorte y alineamiento de miRNA (lo hace nf-core) |
| Matriz de conteos miRNA (isomiR-inclusiva) | piRNA como resultado |
| Exploración (PCA) y DE por contraste (DESeq2) | tRFs y tRNA halves a nivel de fragmento |
| Composición por biotipo (mirTrace; completa desde FASTQ) | miRNAs novel (miRDeep2): módulo opcional de fase 4 |
| Detección de reads derivados de shRNA |  |
| isomiRs, análisis funcional e integración con RNA-seq (fase 3) |  |

Debe servir para otros batches, tejidos y especies: ningún nombre de muestra ni de condición va escrito en el código; todo viene de `metadata.csv`, `contrasts.csv` y parámetros.

## 2. Hechos verificados del batch\_2026\_04

Salen del MultiQC, de las tablas de nf-core y de las comprobaciones hechas sobre los datos. Condicionan el diseño y los tests de aceptación (sección 9) los usan como valores esperados.

| # | Hecho | Consecuencia |
| --- | --- | --- |
| H1 | 49 muestras single-end de 50 nt; 37-90 M reads tras fastp. En crudo, FastQC detecta el adaptador universal de Illumina (AGATCGGAAGAGC…) en 97-99 % de los reads al final de la lectura; tras recortar queda por debajo del 1 %. | El recorte de nf-core funcionó. `run_command.sh` declara el adaptador TruSeq (TGGAATTCTCGGGTGCCAAGG), que no es el presente en los datos (el pico de 22 nt no prueba el adaptador). Los módulos propios usan el adaptador explícito `AGATCGGAAGAGCACACGTCT` (parámetro `adapter_3p`). Confirmar el kit con Novogene. |
| H2 | `mature_counts.csv` (coincidencia exacta con el miRNA canónico) suma 2.141.084 reads en NB4\_Exp7\_Scramble; mirtop suma 22.207.504 (el canónico es el 9,6 %). mirTrace da 20.432.348 reads miRNA. | La cuantificación por defecto es isomiR-inclusiva (mirtop), no `mature_counts`. |
| H3 | En `joined_samples_mirtop.tsv`, 5.202 de 182.547 secuencias únicas (2,8 %) están asignadas a más de un miRNA. `mirna.tsv` agrega sumando, así que cuenta dos veces esos reads (mirtop supera a mirTrace en un 9 %). | Los conteos se calculan desde `joined_samples_mirtop.tsv`, repartiendo cada secuencia 1/k entre sus k miRNAs. `mirna.tsv` solo sirve de cruce. |
| H4 | mirTrace: miRNA 1,5-79 % de los reads (media 44 %), tRNA hasta 52 %, rRNA hasta 33 %, unknown 19-98 %. Longitud mediana de 23 nt en casi todas las muestras y hasta 33 nt en MNC con mucho tRNA. | La composición por biotipo es informativa y variable; los tRNA halves sí se capturan en MNC. |
| H5 | Las 6 muestras DDX41sh1 tienen 1,6-21 % de reads miRNA (NB4: 1,6 / 4,2 / 1,6 %; MSCline: 14,5 / 20,9 / 7,0 %) frente a 53-77 % en sus Scramble. En el primer millón de reads de NB4\_Exp4\_DDX41sh1, \~92 % son variantes (extensiones y recortes en 5' y 3') de TTGATGACGTGCTGGATGTTT. En su Scramble los más abundantes son miRNAs canónicos (miR-148a-3p, miR-21-5p, miR-25-3p, miR-92a-3p…). | Compatible con una guía de shRNA que domina la librería (complemento inverso: AAACATCCAGCACGTCATCAA; comparar con el sitio diana de sh1). Hacen falta M6 (fracción shRNA), flags de composición y análisis de sensibilidad. **\[VERIFICAR\]** con la secuencia real del shRNA. |
| H6 | En `results/` solo hay los `.flagstat`, `.idxstats` y `.stats` de los BAM (no hay BAM), ni FASTQ recortados, ni FastQC; el directorio `work/` se borró. | Los módulos basados en lecturas (M6, M7) parten de los FASTQ crudos de `samplesheet.csv`. El MVP (M0-M5) no necesita relanzar nf-core. |
| H7 | `mature_counts.csv` está traspuesto (muestras en filas, miRNAs en columnas). `mirna.tsv` y `joined_samples_mirtop.tsv` traen las columnas de muestra en orden arbitrario. | Unir siempre por nombre de muestra, nunca por posición. |
| H8 | Los miRNAs de `mature_counts.csv` son humanos (hsa-…, hasta hsa-miR-12136; del orden de 2.600). `mature` y `hairpin` apuntaban a la URL de test-datasets de nf-core. | Es un miRBase completo, no un subconjunto, pero hay que registrar la versión exacta (ver M0). |

## 3. Decisiones de diseño

| ID | Decisión | Motivo |
| --- | --- | --- |
| D1 | Nextflow DSL2 ≥ 24.10 (el clúster usa 25.04.3); plantilla nf-core recomendada (`nf-core pipelines create`), validación de parámetros con nf-schema y tests con nf-test. | Estándar del ecosistema. |
| D2 | Singularity (hay caché compartida en `NXF_SINGULARITY_CACHEDIR`); herramientas de línea de comandos desde biocontainers y una imagen R propia con versiones fijadas. | Reproducibilidad. |
| D3 | `count_source = mirtop` por defecto; `mature` como análisis de sensibilidad. | H2. |
| D4 | Reparto 1/k de las secuencias multi-asignadas; redondeo final a entero. | H3. |
| D5 | Filtros de isomiR provisionales (tabla 5.3) para eliminar asignaciones espurias (SNV en semilla o centro, desplazamientos grandes). **\[VERIFICAR\]** que quitan menos del 10 % de los reads. | Hay secuencias de 17 nt asignadas con 5p +6 o con SNV en la semilla. |
| D6 | Un objeto DESeq2 por `norm_group` (`~ condition`); los contrastes salen de `contrasts.csv` con log2FC = treatment/reference. | Petición del usuario. |
| D7 | `Exp` solo es metadata y color de PCA; no entra en el modelo. | Decisión del usuario. |
| D8 | Las muestras marcadas (flags) no se excluyen; el DE se repite sin ellas como sensibilidad. | Evitar exclusiones silenciosas. |
| D9 | Composición por biotipo en dos niveles: mirTrace (fase 1, inmediato) y completa desde FASTQ con anotación jerárquica (fase 2). | Empezar simple. |
| D10 | M6 y M7 parten de FASTQ crudos con recorte propio y adaptador explícito. | H1, H6. |
| D11 | No se relanza nf-core para el MVP. Si se relanza: misma versión (2.4.1) y mismos parámetros más `--save_intermediates`; no cambiar el adaptador a mitad de batch. | Mantener comparabilidad. |
| D12 | Semillas fijas (`seed = 42`); versiones de software y parámetros volcados en `pipeline_info/`. | Reproducibilidad. |

## 4. Entradas (contratos)

### 4.1 Ficheros de nf-core/smrnaseq que se leen

Rutas relativas a `nfcore_dir` (la carpeta `results/` del run).

| Fichero | Formato | Lo usa | Notas |
| --- | --- | --- | --- |
| `mirna_quant/mirtop/joined_samples_mirtop.tsv` | TSV. Columnas: UID, Read, miRNA, Variant, iso\_5p, iso\_3p, iso\_add3p, iso\_snp, iso\_5p\_nt, iso\_3p\_nt, iso\_add3p\_nt, iso\_snp\_nt, y una columna por muestra (orden arbitrario). Una fila por par (secuencia, miRNA). | M1, M6, M8 | Fuente de conteo por defecto. |
| `mirna_quant/mirtop/mirna.tsv` | TSV miRNA × muestra (sumado; doble conteo en las multi-asignadas) | M1 (solo cruce) | No usar como conteo final. |
| `mirna_quant/edger_qc/mature_counts.csv` | CSV muestras × miRNA (traspuesto); la primera columna, sin nombre, es la muestra | M1 | Se traspone. |
| `mirna_quant/edger_qc/hairpin_counts.csv` | ídem, a nivel de precursor | M1 | Informativo. |
| `multiqc/multiqc_data/multiqc_general_stats.txt` | TSV de MultiQC | M1, M2 | **\[VERIFICAR\]** el nombre exacto con `ls multiqc/multiqc_data`. |
| `mirtrace/<muestra>/` (`mirtrace-results.json`, `mirtrace-stats-*.tsv`) | JSON y TSV de mirTrace | M1, M3 | **\[VERIFICAR\]** claves del JSON y nombres de los TSV (categorías: miRNA, rRNA, tRNA, Artifact, Unknown; distribución de longitudes). |
| `fastp/on_raw/<muestra>.fastp.json` | JSON | M1, M2 | Reads antes y después, longitud. |
| `genome_quant/bam/` (ficheros `.stats`, `.idxstats`, `.flagstat`) | texto de samtools | M2 | Informativo; no hay BAM. |
| `pipeline_info/params_*.json` y `nf_core_smrnaseq_software_mqc_versions.yml` | JSON y YAML | M0 | Auditoría del run. |
| `samplesheet.csv` (carpeta padre de `results/`) | CSV `sample,fastq_1` | M6, M7 | Rutas absolutas a los FASTQ crudos. |

### 4.2 `metadata.csv`

Una fila por muestra.

| Columna | Obligatoria | Descripción |
| --- | --- | --- |
| `sample_id` | sí | Idéntico al nombre de muestra de nf-core (p. ej. NB4\_Exp4\_DDX41sh1). |
| `condition` | sí | Etiqueta que aparece en `contrasts.csv` (p. ej. NB4\_DDX41\_sh1, MNC\_Healthy). |
| `norm_group` | sí | Subconjunto que comparte modelo y normalización (NB4, MSCline, MNC, MSC). |
| `system` | no | NB4, MSCline o primary. |
| `tissue` | no | MNC, MSC o la línea. |
| `exp` | no | Experimento (Exp4…); solo color de PCA. |
| `subject` | no | P1…, H1… |
| `library_id` | no | C1…C49 (identificador del FASTQ). |
| `wts_id` | no | Identificador de la muestra de RNA-seq para la integración (fase 3). |

Generación: `bin/make_metadata.py --samplemap samplemap_rnaseq_2026_04.tsv --samplesheet samplesheet.csv`. Deriva las columnas del nombre de muestra con estas reglas deterministas y contrasta con el samplemap:

- `SISTEMA_ExpN_X`, con SISTEMA igual a NB4 o MSCline y X igual a Scramble, DDX41sh1, DDX41sh2, DHX34sh1 o DHX34sh2: `condition` = `SISTEMA_Scramble` o `SISTEMA_GEN_shK` (p. ej. NB4\_DHX34\_sh2); `norm_group` = SISTEMA.
- `DDX41Patient_TEJIDO_Pn`, con TEJIDO igual a MNC o MSC: `condition` = `TEJIDO_DDX41_Patient`; `norm_group` = TEJIDO.
- `HealthyDonor_TEJIDO_Hn`: `condition` = `TEJIDO_Healthy`; `norm_group` = TEJIDO.

Un nombre que no case con ninguna regla detiene el script (no se adivina).

### 4.3 `contrasts.csv` (formato fijado)

```
contrast_id,variable,reference,treatment,norm_group
MSCline_DDX41_sh1_vs_Scramble,condition,MSCline_Scramble,MSCline_DDX41_sh1,MSCline
MSCline_DDX41_sh2_vs_Scramble,condition,MSCline_Scramble,MSCline_DDX41_sh2,MSCline
MSCline_DHX34_sh1_vs_Scramble,condition,MSCline_Scramble,MSCline_DHX34_sh1,MSCline
MSCline_DHX34_sh2_vs_Scramble,condition,MSCline_Scramble,MSCline_DHX34_sh2,MSCline
NB4_DDX41_sh1_vs_Scramble,condition,NB4_Scramble,NB4_DDX41_sh1,NB4
NB4_DDX41_sh2_vs_Scramble,condition,NB4_Scramble,NB4_DDX41_sh2,NB4
NB4_DHX34_sh1_vs_Scramble,condition,NB4_Scramble,NB4_DHX34_sh1,NB4
NB4_DHX34_sh2_vs_Scramble,condition,NB4_Scramble,NB4_DHX34_sh2,NB4
MNC_DDX41_Patient_vs_Healthy,condition,MNC_Healthy,MNC_DDX41_Patient,MNC
MSC_DDX41_Patient_vs_Healthy,condition,MSC_Healthy,MSC_DDX41_Patient,MSC
```

Validación (M0): `variable` es una columna de `metadata.csv`; `reference` y `treatment` existen en esa columna y pertenecen al `norm_group` indicado; `contrast_id` es único; `reference` y `treatment` difieren; hay al menos 2 muestras por nivel (con menos de 3, aviso). Signo: log2FC positivo significa más abundante en `treatment` que en `reference`; log2FC = 2 equivale a ×4.

### 4.4 Referencias (por proyecto)

| Recurso | Para qué | Estado |
| --- | --- | --- |
| FASTA del genoma (GRCh38 Ensembl, `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`) | M7 | existe |
| Índice Bowtie1 del genoma | M7 | **\[VERIFICAR\]** `data/reference/genome/GRCh38/bowtie1_index/`; si no existe, se construye desde el FASTA |
| GTF de Ensembl (misma release que el FASTA) | M7 | pendiente de indicar |
| miRBase hsa.gff3 (GRCh38) | M7 | pendiente |
| Anotación de tRNA nuclear (tRNAs de GENCODE o gtRNAdb hg38) | M7 | pendiente; hasta donde sé, el GTF de Ensembl solo trae los tRNA mitocondriales |
| miRTarBase (fichero local) | M9 | opcional |

Nombres de cromosoma: Ensembl usa 1 y MT; GENCODE y miRBase usan chr1 y chrM. M7 incluye una tabla de equivalencia (`assets/chr_map.tsv`) y falla si algún contig no se resuelve.

## 5. Parámetros

Tres niveles. **A**: cambian en cada ejecución (van en el `params.yaml` del run). **B**: cambian por proyecto, kit o especie (van en un `params.yaml` o config del proyecto). **C**: defaults del pipeline (en `nextflow.config`); solo se tocan por una decisión de análisis.

### 5.1 Nivel A: por ejecución

| Parámetro | Default | Descripción |
| --- | --- | --- |
| `nfcore_dir` | (obligatorio) | Carpeta `results/` de nf-core/smrnaseq. |
| `metadata` | (obligatorio) | `metadata.csv`. |
| `contrasts` | (obligatorio) | `contrasts.csv`. |
| `outdir` | `results_smrna` | Carpeta de salida. |
| `fastq_samplesheet` | `nfcore_dir/../samplesheet.csv` | Solo si se activan M6 o M7. |
| `shrna_fasta` | null | FASTA con la secuencia guía de cada shRNA (cabecera = nombre). Opcional pero muy recomendable. |
| `consistency_pairs` | null | CSV `pair_id,contrast_a,contrast_b` con los contrastes cuya concordancia se compara (p. ej. sh1 frente a sh2). |
| `rnaseq_de_dir`, `integration_map` | null | Resultados de DE de RNA-seq y tabla de equivalencia de contrastes (fase 3). |
| `run_shrna`, `run_biotype_full`, `run_isomir`, `run_functional` | false | Activan los módulos de las fases 2 y 3. M0-M5 y M10 siempre se ejecutan. |

### 5.2 Nivel B: por proyecto, kit o especie

| Parámetro | Default | Descripción |
| --- | --- | --- |
| `species` | `hsa` | Prefijo de miRNA y especie de mirTrace. |
| `genome_build` | `GRCh38` | Etiqueta; se registra en la auditoría. |
| `genome_fasta`, `genome_bowtie_index`, `gtf`, `mirbase_gff3`, `trna_annotation`, `chr_map` | null | Referencias de la sección 4.4 (M7). |
| `adapter_3p` | `AGATCGGAAGAGCACACGTCT` | Adaptador 3' del kit (M6 y M7). |
| `mirna_db` | `mirbase` | `mirbase` o `mirgenedb`; cambia la armonización de nombres y desactiva M9 con MirGeneDB. |
| `mirtarbase_file` | null | miRTarBase local (M9). |

### 5.3 Nivel C: defaults del pipeline

| Parámetro | Default | Descripción |
| --- | --- | --- |
| `count_source` | `mirtop` | `mirtop` o `mature`. |
| `multi_assign` | `split` | `split` (1/k), `all` o `drop`. |
| `round_counts` | true | Redondear a entero tras sumar. |
| `isomir_exclude_variants` | `iso_snv_seed,iso_snv_central,iso_snv_central_offset,iso_snv_central_supp` | Las filas con alguna de estas etiquetas en `Variant` se descartan. |
| `isomir_max_shift_5p`, `isomir_max_shift_3p` | 2, 5 | Valor absoluto máximo de iso\_5p e iso\_3p. |
| `isomir_min_len` | 18 | Longitud mínima del read. |
| `qc_min_mirna_frac` | 0.10 | Flag LOW\_MIRNA\_FRAC si la fracción miRNA (mirTrace) es menor. |
| `qc_max_trna_frac`, `qc_max_rrna_frac` | 0.30, 0.20 | Flags HIGH\_TRNA\_FRAC y HIGH\_RRNA\_FRAC. |
| `qc_min_mirna_reads` | 1000000 | Flag LOW\_DEPTH (suma de conteos mirtop). |
| `sensitivity_exclude_flags` | `LOW_MIRNA_FRAC,LOW_DEPTH` | Flags que definen las muestras excluidas en el DE de sensibilidad. |
| `min_cpm`, `min_samples`, `min_total_count` | 1, auto, 10 | Filtro de expresión: CPM ≥ `min_cpm` en ≥ `min_samples` muestras (auto = tamaño del nivel más pequeño del grupo) y suma ≥ `min_total_count`. |
| `alpha`, `lfc_flag` | 0.05, 1 | FDR y umbral de valor absoluto de log2FC para marcar efecto fuerte. |
| `shrink_type` | `ashr` | Contracción de log2FC (`lfcShrink` con contraste). |
| `sizefactor_type` | `ratio` | Con paso automático a `poscounts` si falla (se registra en el log). |
| `n_top_var`, `n_top_heatmap` | 500, 30 | miRNAs para la PCA y para el heatmap de DE. |
| `sensitivity_runs` | `mature,no_flagged` | Repeticiones del DE. |
| `dominant_seq_frac` | 0.20 | Aviso si una secuencia o familia ocupa más de esta fracción de los reads (M6). |
| `top_n_sequences` | 50 | Secuencias que lista M6 por muestra. |
| `shrna_prefix_len`, `shrna_max_mismatch` | 19, 1 | Criterio de coincidencia con la guía (longitud del prefijo y desajustes). |
| `trim_min_len`, `trim_max_len` | 17, 100 | Longitudes tras fastp (igual que nf-core). |
| `bowtie_mismatches`, `bowtie_max_hits` | 1, 50 | M7: desajustes permitidos y límite de alineamientos por secuencia. |
| `min_overlap_frac` | 0.8 | Fracción mínima del read dentro de una feature para anotarlo. |
| `biotype_priority` | `miRNA,tRNA,rRNA,snoRNA,snRNA,scaRNA,Y_RNA,vault_RNA,other_ncRNA,lncRNA,protein_coding,unannotated` | Orden de la asignación jerárquica. |
| `de_classes` | `snoRNA,snRNA,tRNA,rRNA` | Clases con DE a nivel de feature (M7). |
| `exclude_shrna_from_composition` | true | Con `shrna_fasta`, la composición se calcula sin esos reads (y también con ellos). |
| `isomir_min_reads` | 500 | Lecturas mínimas por miRNA y muestra para el perfil de isomiRs (M8). |
| `min_sig_mirnas` | 5 | Mínimo de miRNAs significativos para ejecutar el enriquecimiento (M9). |
| `seed` | 42 | Semilla. |
| `plot_formats`, `plot_dpi` | `pdf,png`, 300 | Formatos de figura. |

### 5.4 `params.yaml` de ejemplo (este batch)

```
nfcore_dir: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/smRNAseq/batch_2026_04/nfcore_smrnaseq/results
metadata: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/smRNAseq/batch_2026_04/smrna_nextflow/metadata.csv
contrasts: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/smRNAseq/batch_2026_04/smrna_nextflow/contrasts.csv
consistency_pairs: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/smRNAseq/batch_2026_04/smrna_nextflow/consistency_pairs.csv
outdir: results_smrna
species: hsa
genome_fasta: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/data/reference/genome/GRCh38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
adapter_3p: AGATCGGAAGAGCACACGTCT
run_shrna: false
run_biotype_full: false
```

### 5.5 Lanzador y perfiles

- Código en `CRUK_DDX41_DHX34/pipelines/small_RNAseq/smrna_nextflow/` (junto a `nf-core_smrnaseq/`); ejecución en `analysis/smRNAseq/BATCH/smrna_nextflow/`.
- Perfiles: `singularity`, `slurm` (executor slurm, cola `compute`) y `test`. La config institucional va en `conf/custom.config` (misma lógica que la del run de nf-core).
- `run_smrna.sh` replica el patrón de `run_command.sh`: cargar el módulo de nextflow, fijar `NXF_VER` y `NXF_SINGULARITY_CACHEDIR`, `-work-dir` en scratch (nunca dentro del proyecto) y `-resume`. El job de cabecera necesita poco (2 CPU, 8 GB), igual que `submit.sh`.

## 6. Arquitectura

```
smrna_nextflow/
├── main.nf
├── nextflow.config, nextflow_schema.json
├── conf/{base,modules,slurm,test,custom}.config
├── workflows/smrna_nextflow.nf
├── subworkflows/local/{ingest_qc,differential,reads_based}/main.nf
├── modules/local/MODULO/{main.nf,meta.yml,environment.yml}
├── bin/            # scripts R y Python (ejecutables, con --help)
├── assets/         # report.qmd, chr_map.tsv, biotype_map.tsv, esquemas
├── docs/           # README, uso, esta especificación
├── tests/          # nf-test y datos de prueba (6 muestras)
└── Dockerfile, renv.lock   # imagen R
```

Flujo por fases:

```
Fase 1 (MVP):  M0 VALIDATE -> M1 INGEST -> M2 QC -> M3 BIOTYPE_MIRTRACE
                                      \-> M4 EXPLORE -> M5 DE -> M10 REPORT
Fase 2:        TRIM_COLLAPSE -> M6 SHRNA_TOPSEQ
                             \-> M7 BIOTYPE_FULL
Fase 3:        M8 ISOMIR, M9 FUNCTIONAL
Fase 4 (opc.): M11 NOVEL_MIRNA
```

Contenedores: biocontainers para fastp, vsearch, bowtie, samtools, bedtools y seqkit; una imagen R (R ≥ 4.4 con la versión de Bioconductor correspondiente) con DESeq2, ashr, ggplot2, ComplexHeatmap, data.table, optparse, jsonlite, yaml y quarto, y en la fase 3 clusterProfiler y fgsea. Las versiones exactas se fijan en `renv.lock` en el primer build y se vuelcan a `pipeline_info/software_versions.yml`. Etiquetas de recursos: `process_single` (M0, M1, M2, M3), `process_medium` (M4, M5, M8) y `process_high` (TRIM\_COLLAPSE y M7: 8 CPU y 32 GB como punto de partida).

## 7. Especificación de módulos

Cada módulo indica entradas, pasos, salidas con su contrato de columnas y condiciones de fallo. Los umbrales citados son los defaults de la tabla 5.3.

### M0 · VALIDATE\_INPUTS (fase 1)

- **Entradas:** `nfcore_dir`, `metadata`, `contrasts`, `pipeline_info/params_*.json`, `pipeline_info/nf_core_smrnaseq_software_mqc_versions.yml`.
- **Pasos:**
  1. Comprobar que existen todos los ficheros de la tabla 4.1 (con glob) y listarlos en `inputs_manifest.tsv` (ruta, tamaño; sha256 solo de metadata y contrasts).
  2. Leer `params_*.json` y registrar: versión del pipeline, genoma, `mature` y `hairpin`, `mirgenedb`, adaptador declarado, `fastp_min_length`, `with_umi`, `filter_contamination`. Emitir avisos (no fallos) si: el adaptador declarado difiere de `adapter_3p`; `mirgenedb` no coincide con `mirna_db`; `with_umi` es true; `filter_contamination` es true (cambia la lectura de las composiciones).
  3. Validar `metadata.csv`: columnas obligatorias, `sample_id` únicos, todas las muestras presentes en `mature_counts.csv` y en `joined_samples_mirtop.tsv`, ninguna muestra sobrante.
  4. Validar `contrasts.csv` (sección 4.3).
  5. Registrar el número de muestras por condición y por `norm_group`.
- **Salidas** (`00_audit/`): `audit.json`, `audit.md` (legible), `inputs_manifest.tsv`, `metadata_validated.csv`, `contrasts_validated.csv`.
- **Fallo duro:** ficheros obligatorios ausentes, muestras que no casan, contraste inválido. El resto son avisos que pasan al informe.

### M1 · INGEST\_COUNTS (fase 1)

- **Entradas:** los ficheros de la tabla 4.1 y `metadata_validated.csv`.
- **Pasos (fuente mirtop):**
  1. Leer `joined_samples_mirtop.tsv` (`data.table::fread`) y quedarse con las columnas de muestra de `metadata`.
  2. Aplicar los filtros de isomiR (D5): descartar filas cuyo `Variant` contenga alguna etiqueta de `isomir_exclude_variants`, filas con valor absoluto de `iso_5p` mayor que `isomir_max_shift_5p` o de `iso_3p` mayor que `isomir_max_shift_3p` (los valores vienen como +1, -2, 0: convertir a entero) y reads con longitud menor que `isomir_min_len`. Contar los reads retirados por muestra.
  3. Calcular k = número de miRNAs distintos por UID tras el filtrado. Peso = 1/k con `split`; 1 con `all`; 0 si k > 1 con `drop`.
  4. Conteo(muestra, miRNA) = suma de peso × conteo. Guardar la versión fraccionaria y redondear a entero (`round_counts`).
- **Pasos (fuente mature):** leer `mature_counts.csv`, trasponer y ordenar las columnas según `metadata`.
- **Pasos comunes:** 5. Armonizar nombres de miRNA entre fuentes (coincidencia exacta); los no emparejados van a `mirna_id_mismatch.tsv` (no es un error). 6. Cruce de totales por muestra: joined sin filtrar, `mirna.tsv`, mirtop final, mature y mirTrace. 7. Construir `qc_table.tsv` fusionando fastp (reads antes y después, % que pasan, longitud media), mirTrace (reads y % por categoría), samtools (% mapeado), totales de mirtop y metadata.
- **Salidas:** `01_counts/counts_mirtop.tsv` (miRNA × muestra, enteros, columnas en el orden de `metadata`), `counts_mirtop_fractional.tsv`, `counts_mature.tsv`, `counts_hairpin.tsv`, `isomir_filtered.tsv.gz` (filas supervivientes con su peso), `ingest_stats.tsv`, `mirna_id_mismatch.tsv`; y `02_qc/qc_table.tsv`.
- **Columnas de `ingest_stats.tsv`:** `sample_id`, `joined_raw_total`, `mirna_tsv_total`, `rows_removed_frac`, `mirtop_final_total`, `mature_total`, `mirtrace_mirna_reads`, `canonical_frac` (= `mature_total` / `mirtop_final_total`), `multi_assigned_frac`.
- **Columnas de `qc_table.tsv`:** `sample_id`, `condition`, `norm_group`, `exp`, `reads_raw`, `reads_passed_fastp`, `pct_passed`, `mean_len_trimmed`, `mirtrace_pct_miRNA`, `mirtrace_pct_tRNA`, `mirtrace_pct_rRNA`, `mirtrace_pct_artifact`, `mirtrace_pct_unknown`, `mirna_reads_mirtop`, `canonical_frac`, `pct_dup`, `genome_mapped_pct`.
- **Validaciones:** conteos enteros no negativos y sin NA; orden de columnas igual al de `metadata`; `mirtop_final_total` no mayor que `joined_raw_total`.

### M2 · QC\_SUMMARY (fase 1)

- **Entrada:** `qc_table.tsv`.
- **Flags** (se marcan, no se excluyen): LOW\_MIRNA\_FRAC (fracción miRNA de mirTrace menor que `qc_min_mirna_frac`), HIGH\_TRNA\_FRAC, HIGH\_RRNA\_FRAC, LOW\_DEPTH (`mirna_reads_mirtop` menor que `qc_min_mirna_reads`) y, desde la fase 2, DOMINANT\_SEQUENCE (M6).
- **Salidas** (`02_qc/`): `qc_flags.tsv` (`sample_id`, `flag`, `value`, `threshold`), `flagged_samples.txt` (muestras con algún flag de `sensitivity_exclude_flags`), `qc_overview.pdf` y `.png` (reads por muestra, % miRNA, longitud mediana y fracción canónica, coloreados por condición).

### M3 · BIOTYPE\_MIRTRACE (fase 1)

- **Entradas:** `qc_table.tsv`, metadata y contrasts.
- **Pasos:** porcentaje por muestra de cada categoría de mirTrace (miRNA, rRNA, tRNA, Artifact, Unknown). Por contraste: media por condición y Δ = media(treatment) − media(reference) en puntos porcentuales. Prueba indicativa: transformación CLR (pseudoconteo de 0,5 reads) y t de Welch por categoría con FDR de Benjamini-Hochberg entre categorías; `indicative_only = TRUE` si algún grupo tiene menos de 4 muestras.
- **Salidas** (`03_biotype_mirtrace/`): `composition_by_sample.tsv`, `delta_by_contrast.tsv` (`contrast_id`, `category`, `mean_ref`, `mean_trt`, `delta_pp`, `t_stat`, `p`, `fdr`, `indicative_only`), `composition_stacked_bar.pdf` (muestras ordenadas por `norm_group` y `condition`) y `delta_CONTRAST_ID.pdf`.
- **Nota:** los reads derivados de shRNA caen en Unknown hasta que M6 y M7 los separan; el informe lo avisa.

### M4 · FILTER\_NORM\_EXPLORE (fase 1)

- **Entradas:** conteos (según `count_source`), metadata y `qc_flags.tsv`.
- **Pasos** (por `norm_group`, más un análisis global de las 49 muestras sin DE):
  1. Subconjunto de muestras del grupo.
  2. Filtro de expresión (`min_cpm`, `min_samples`, `min_total_count`).
  3. DESeqDataSet `~ condition`; `estimateSizeFactors` con `sizefactor_type` y paso automático a `poscounts` si falla.
  4. VST con `blind = TRUE` para la exploración.
  5. PCA con los `n_top_var` miRNAs de mayor varianza: PC1-PC2 y PC2-PC3, coloreada por `condition` con forma por `exp`; una segunda figura coloreada por `system` o `tissue` y una tercera por flags de QC.
  6. Mapa de calor de distancias euclídeas entre muestras con anotación de `condition`, `exp` y flags.
  7. Size factors y su correlación con la fracción miRNA.
- **Salidas** (`04_explore/GRUPO/` y `04_explore/global/`): `pca_condition.pdf`, `pca_exp.pdf`, `pca_flags.pdf`, `pca_coords.tsv`, `sample_distance_heatmap.pdf`, `size_factors.tsv`, `vst_counts.tsv`, `filtered_counts.tsv`, `filter_stats.tsv` (`group`, `n_features_input`, `n_features_kept`, `min_samples_used`).

### M5 · DE\_DESEQ2 (fase 1)

- **Entradas:** conteos filtrados por grupo, contrasts, metadata y `flagged_samples.txt`.
- **Ejecuciones:** `main` (fuente `count_source`, todas las muestras), `sens_mature` (fuente mature) y `sens_no_flagged` (sin muestras marcadas; un contraste con menos de 2 muestras en un nivel se omite y se registra).
- **Pasos:**
  1. Por `norm_group`: DESeqDataSet `~ condition`.
  2. `DESeq(dds)` con `fitType` paramétrico; si falla, `local` y luego `mean` (se registra).
  3. Por contraste: `results(dds, contrast = c(variable, treatment, reference), alpha = alpha)` y `lfcShrink(dds, contrast = ..., type = shrink_type)`.
  4. Clasificar: `direction` = up (padj menor que `alpha` y log2FC positivo), down o ns; `strong` si el valor absoluto de log2FC es al menos `lfc_flag`.
  5. Concordancia entre `main` y cada sensibilidad: Spearman y Pearson de log2FC (MLE) sobre los miRNAs comunes y número de miRNAs significativos compartidos.
  6. Con `consistency_pairs`: por cada pareja, concordancia de signo y dispersión de log2FC entre ambos contrastes (comprobación sh1 frente a sh2).
- **Salidas por ejecución y contraste** (`05_de/RUN/CONTRAST_ID/`): `de_results.tsv` con columnas `mirna`, `baseMean`, `log2FC_mle`, `log2FC_shrunk`, `lfcSE`, `stat`, `pvalue`, `padj`, `direction`, `strong`, `n_ref`, `n_trt`, `mean_norm_count_ref`, `mean_norm_count_trt`; `ma_plot.pdf`; `volcano.pdf` (etiquetas de los 10 mejores por padj); `heatmap_top.pdf` (los `n_top_heatmap` mejores por padj, z-score de VST por fila).
- **Resúmenes** (`05_de/`): `de_summary.tsv` (`run`, `contrast_id`, `n_tested`, `n_up`, `n_down`, `n_strong`), `concordance.tsv` (`contrast_id`, `run_a`, `run_b`, `n_common`, `spearman_log2fc`, `n_sig_shared`), `concordance_CONTRAST_ID.pdf` y `shRNA_consistency.tsv` (`pair_id`, `n_common`, `spearman_log2fc`, `n_sig_either`, `n_same_sign`).
- **Fichero propuesto para este batch** (`consistency_pairs.csv`):

```
pair_id,contrast_a,contrast_b
NB4_DDX41,NB4_DDX41_sh1_vs_Scramble,NB4_DDX41_sh2_vs_Scramble
NB4_DHX34,NB4_DHX34_sh1_vs_Scramble,NB4_DHX34_sh2_vs_Scramble
MSCline_DDX41,MSCline_DDX41_sh1_vs_Scramble,MSCline_DDX41_sh2_vs_Scramble
MSCline_DHX34,MSCline_DHX34_sh1_vs_Scramble,MSCline_DHX34_sh2_vs_Scramble
```

### TRIM\_COLLAPSE (fase 2; proceso compartido por M6 y M7)

- **Entrada:** FASTQ crudo de cada muestra (`samplesheet.csv`).
- **Pasos:**
  1. `fastp` en modo single-end con `--adapter_sequence` igual a `adapter_3p`, `--length_required` igual a `trim_min_len`, `--length_limit` igual a `trim_max_len` y el resto de valores por defecto de fastp (comparable con nf-core).
  2. Colapso a secuencias únicas con `vsearch --fastx_uniques --sizeout --relabel seq --minseqlength 17`; salida FASTA con cabecera `seqN;size=CONTEO`.
- **Salidas** (`07_reads/`): `MUESTRA.collapsed.fa.gz` y `trim_stats.tsv` (`sample_id`, `reads_in`, `reads_pass`, `unique_seqs`). Intermedio borrable.
- **Recursos:** `process_high`.

### M6 · SHRNA\_TOPSEQ (fase 2, `run_shrna`)

- **Entradas:** `*.collapsed.fa.gz`, `shrna_fasta` (opcional) y `joined_samples_mirtop.tsv` (columna Read, solo como lista de secuencias de miRNA conocidas).
- **Pasos:**
  1. Por muestra, tomar las `top_n_sequences` secuencias con más reads y su porcentaje sobre el total.
  2. Clasificar cada una como `miRNA` (la secuencia figura en la columna Read de mirtop), `shRNA:NOMBRE` (algún tramo de `shrna_prefix_len` nt que empiece en las posiciones 0 a 2 del read coincide con los primeros nt de la guía con a lo sumo `shrna_max_mismatch` desajustes), `shRNA_passenger:NOMBRE` (igual, pero frente al complemento inverso de la guía) o `unclassified`.
  3. Fracción shRNA por muestra = reads de las secuencias que coinciden con cada guía / reads totales.
  4. Familia dominante (con o sin `shrna_fasta`): agrupar las secuencias top que comparten un tramo de 19 nt con desplazamiento de hasta 2 nt; si una familia supera `dominant_seq_frac`, flag DOMINANT\_SEQUENCE.
- **Salidas** (`06_shrna/`): `top_sequences/MUESTRA.top50.tsv` (`rank`, `sequence`, `count`, `pct_reads`, `class`, `matched_name`), `shrna_fraction.tsv` (`sample_id`, `reads_total`, una columna `reads_NOMBRE` por guía, `frac_shrna_total`, `dominant_family`, `dominant_frac`, `flag`) y `dominant_sequences.pdf` (barras apiladas shRNA, miRNA y otros por muestra). Los flags se añaden a `qc_flags.tsv`.
- **Resultado esperado en este batch:** la familia TTGATGACGTGCTGGATGTTT supera 0,85 de los reads en NB4\_Exp4\_DDX41sh1 y ninguna familia supera 0,05 en los Scramble. **\[VERIFICAR\]** con las secuencias reales de los shRNA.

### M7 · BIOTYPE\_FULL (fase 2, `run_biotype_full`)

- **Entradas:** `*.collapsed.fa.gz`, índice Bowtie1 (o `genome_fasta`), `gtf`, `mirbase_gff3`, `trna_annotation`, `chr_map`, metadata, contrasts y, opcionalmente, `shrna_fraction.tsv`.
- **Pasos:**
  1. Preparar la anotación (una vez): `annotation.bed` (`chr`, `start`, `end`, `feature_id`, `class`, `strand`) a partir del GTF, miRBase y tRNA, con nombres de cromosoma armonizados. El paso de `gene_biotype` a clase usa `assets/biotype_map.tsv` **\[OBLIGATORIO\]**: miRNA a miRNA; Mt\_tRNA y tRNA a tRNA; rRNA, Mt\_rRNA y rRNA\_pseudogene a rRNA; snoRNA a snoRNA; snRNA a snRNA; scaRNA a scaRNA; vault\_RNA a vault\_RNA; misc\_RNA cuyo nombre empiece por Y\_RNA o RNY a Y\_RNA; otros ncRNA pequeños a other\_ncRNA; lncRNA y similares (lincRNA, antisense, processed\_transcript, sense\_intronic, sense\_overlapping) a lncRNA; protein\_coding a protein\_coding. **\[VERIFICAR\]** los biotipos reales del GTF con `grep -o 'gene_biotype .[^;]*' GTF | sort | uniq -c` antes de cerrar la tabla.
  2. Alinear con `bowtie -f -v bowtie_mismatches -k (bowtie_max_hits + 1) --best --strata -S`. Las secuencias sin alineamiento van a la clase `unaligned`; las que superan `bowtie_max_hits` alineamientos, a `multimap_gt_max`.
  3. Pasar a BED conservando el nombre y el tamaño de la secuencia y anotar con `bedtools intersect -wo -s -f` igual a `min_overlap_frac` (misma hebra: las librerías de small RNA son sentido). Un alineamiento sin feature solapada va a `unannotated`.
  4. Asignar la clase de mayor prioridad según `biotype_priority` entre las features solapadas; dentro de una clase, la de mayor solapamiento (empate: `feature_id` por orden alfabético).
  5. Ponderar: una secuencia de tamaño n con N loci aporta n/N a cada locus; la suma por clase da `biotype_counts`.
  6. Con `shrna_fraction.tsv` y `exclude_shrna_from_composition`, repetir la composición sin los reads que coinciden con las guías.
  7. Para cada clase de `de_classes`, matriz feature × muestra (tRNA agrupado por anticodón: los tres primeros campos del nombre, p. ej. tRNA-Gly-GCC **\[VERIFICAR\]**) y DE con el mismo código que M5.
- **Salidas** (`07_biotype_full/`): `alignment_stats.tsv`, `biotype_counts.tsv` (clase × muestra, fraccionario), `biotype_composition.tsv` (porcentaje), `delta_by_contrast.tsv` (mismas columnas que M3), `composition_stacked_bar.pdf`, `delta_CONTRAST_ID.pdf`, `feature_counts/CLASE.tsv`, `de/CLASE/CONTRAST_ID/de_results.tsv` y las variantes sin shRNA (sufijo `_noshrna`).
- **Advertencia para el informe:** los porcentajes son composicionales: si una clase aumenta, las demás bajan. Interpretar siempre junto a los conteos absolutos; los DE de clases con pocas features (p. ej. snRNA) son exploratorios.
- **Recursos:** `process_high`.

### M8 · ISOMIR (fase 3, `run_isomir`)

- **Entradas:** `isomir_filtered.tsv.gz`, metadata y contrasts.
- **Pasos:**
  1. Por muestra y miRNA con al menos `isomir_min_reads` lecturas: `reads_total`, `reads_canonical` (iso\_5p igual a 0, iso\_3p igual a 0, sin adición 3' ni SNV), `frac_canonical`, `frac_5p_shifted`, `frac_3p_shifted`, `frac_add3p` e isoforma dominante.
  2. Matriz de isomiRs: feature = miRNA más variante (iso\_5p, iso\_3p, iso\_add3p) unidos con dos puntos, p. ej. `hsa-miR-21-5p:0:+2:CA`; filtro `min_cpm` y DE como en M5.
  3. Por contraste: cambio de `frac_5p_shifted` por miRNA (t de Welch sobre logit, indicativo). Un cambio de desplazamiento 5' altera la semilla y se destaca.
- **Salidas** (`08_isomir/`): `isomir_profile_by_sample.tsv`, `isomir_counts.tsv`, `de_isomir/CONTRAST_ID/de_results.tsv`, `delta_5p_shift_CONTRAST_ID.tsv` y figuras (mapa de calor de fracciones por miRNA; barras de isoformas de los miRNAs más abundantes).

### M9 · FUNCTIONAL (fase 3, `run_functional`)

- **Entradas:** `de_results.tsv` (ejecución `main`), `mirtarbase_file` y, opcionalmente, `rnaseq_de_dir` e `integration_map`.
- **Pasos:**
  1. Targets validados: filtrar miRTarBase a `species` y evidencia funcional; casar por nombre de miRNA (misma versión de miRBase) y reportar los no emparejados.
  2. Por contraste y dirección (up y down por separado): si hay al menos `min_sig_mirnas` miRNAs significativos, enriquecimiento (GO Biological Process) con clusterProfiler sobre sus targets; el universo se define explícitamente (genes con algún target validado de un miRNA expresado).
  3. Integración (solo con `rnaseq_de_dir`): `integration_map.csv` con `contrast_smrna`, `contrast_rnaseq` y `de_file`; el DE de RNA-seq debe ser un TSV con `gene_symbol`, `log2FoldChange`, `stat` y `padj`. Para cada miRNA significativo con al menos 10 targets validados, fgsea de sus targets sobre el ranking `stat`; si el miRNA sube, se espera NES negativo.
- **Salidas** (`09_functional/CONTRAST_ID/`): `targets.tsv`, `ora_up.tsv`, `ora_down.tsv`, `integration_fgsea.tsv` e `integration_plot.pdf`.
- Fuera de alcance: predicción de targets por secuencia.

### M10 · REPORT (fase 1, se amplía por fases)

- Informe Quarto en HTML que solo lee `results_smrna/` (no recalcula nada) y falla si falta un fichero obligatorio de las fases activas.
- **Secciones:** 1 resumen ejecutivo (muestras, flags, contrastes y nº de miRNAs diferenciales); 2 auditoría del run; 3 QC y flags; 4 composición (mirTrace y, si existe, completa); 5 shRNA y secuencias dominantes (si existe); 6 exploración (PCA); 7 DE por contraste con la concordancia entre main, mature y sin muestras marcadas; 8 isomiRs y análisis funcional (si existen); 9 métodos (herramientas, versiones, parámetros); 10 limitaciones (texto fijo: composición afectada por reads de shRNA; el conteo canónico solo recoge una fracción de los reads; n pequeño).
- **Salidas** (`report/`): `smrna_report.html`, `tables/` (copias en TSV de las tablas principales) y `figures/`.

### M11 · NOVEL\_MIRNA (fase 4, opcional)

No implementar hasta cerrar las fases 1-3. Esquema: leer `mirdeep2/result_MUESTRA.csv` y `.bed` **\[VERIFICAR el formato\]**; quedarse con los candidatos novel por encima de un score mínimo de miRDeep2 y con lectura madura suficiente; fusionar por solapamiento de coordenadas entre muestras (en al menos 3); descartar solapes con miRBase, tRNA y snoRNA; recuantificar. Salida: `11_novel/novel_candidates.tsv`, siempre etiquetada como candidata.

## 8. Estructura de outputs y entregables

```
results_smrna/
├── 00_audit/             audit.json, audit.md, inputs_manifest.tsv, metadata_validated.csv, contrasts_validated.csv
├── 01_counts/            counts_mirtop.tsv, counts_mirtop_fractional.tsv, counts_mature.tsv, counts_hairpin.tsv,
│                         isomir_filtered.tsv.gz, ingest_stats.tsv, mirna_id_mismatch.tsv
├── 02_qc/                qc_table.tsv, qc_flags.tsv, flagged_samples.txt, qc_overview.pdf y .png
├── 03_biotype_mirtrace/  composition_by_sample.tsv, delta_by_contrast.tsv, composition_stacked_bar.pdf, delta_CONTRAST_ID.pdf
├── 04_explore/           global/ y GRUPO/: pca_*.pdf, pca_coords.tsv, sample_distance_heatmap.pdf,
│                         size_factors.tsv, vst_counts.tsv, filtered_counts.tsv, filter_stats.tsv
├── 05_de/                main/CONTRAST_ID/{de_results.tsv, ma_plot.pdf, volcano.pdf, heatmap_top.pdf}
│                         sens_mature/..., sens_no_flagged/..., de_summary.tsv, concordance.tsv,
│                         concordance_CONTRAST_ID.pdf, shRNA_consistency.tsv (si hay consistency_pairs)
├── 06_shrna/             top_sequences/MUESTRA.top50.tsv, shrna_fraction.tsv, dominant_sequences.pdf     [fase 2]
├── 07_reads/             MUESTRA.collapsed.fa.gz, trim_stats.tsv (intermedio, borrable)                  [fase 2]
├── 07_biotype_full/      alignment_stats.tsv, biotype_counts.tsv, biotype_composition.tsv,
│                         delta_by_contrast.tsv, feature_counts/, de/CLASE/CONTRAST_ID/, figuras         [fase 2]
├── 08_isomir/            isomir_profile_by_sample.tsv, isomir_counts.tsv, de_isomir/, figuras             [fase 3]
├── 09_functional/        CONTRAST_ID/{targets.tsv, ora_up.tsv, ora_down.tsv, integration_fgsea.tsv}      [fase 3]
├── report/               smrna_report.html, tables/, figures/
└── pipeline_info/        params.json, software_versions.yml, execution_report.html, execution_timeline.html,
                          execution_trace.txt, pipeline_dag.html
```

### Entregables: qué abrir para cada pregunta

| Pregunta | Archivo |
| --- | --- |
| ¿Están bien las muestras? | `02_qc/qc_overview.pdf`, `02_qc/qc_flags.tsv` |
| ¿Cómo es la composición del small RNA? | `03_biotype_mirtrace/composition_stacked_bar.pdf` (y, en fase 2, `07_biotype_full/composition_stacked_bar.pdf`) |
| ¿Cómo cambia la composición en cada contraste? | `03_biotype_mirtrace/delta_by_contrast.tsv` y `delta_CONTRAST_ID.pdf` |
| ¿Se agrupan las muestras por biología? | `04_explore/GRUPO/pca_condition.pdf` |
| ¿Qué miRNAs cambian en cada contraste? | `05_de/main/CONTRAST_ID/de_results.tsv` y `volcano.pdf` |
| ¿Es robusto el resultado? | `05_de/concordance.tsv` y `05_de/shRNA_consistency.tsv` |
| ¿Hay una secuencia (shRNA) dominando la librería? | `06_shrna/shrna_fraction.tsv` |
| ¿Cambia el procesamiento de los miRNAs? | `08_isomir/delta_5p_shift_CONTRAST_ID.tsv` |
| Resumen para compartir | `report/smrna_report.html` |

Convenciones de ficheros: snake\_case; TSV con cabecera, separador tabulador, sin comillas y con punto decimal; `sample_id` o `contrast_id` como primera columna cuando aplique; miRNAs con nombre completo (hsa-miR-21-5p); figuras en PDF (vectorial) y PNG.

## 9. Criterios de aceptación y tests

### 9.1 Valores esperados en este batch

| ID | Test | Esperado |
| --- | --- | --- |
| A1 | Suma de `counts_mature.tsv` en NB4\_Exp7\_Scramble | 2.141.084 |
| A2 | `joined_raw_total` de NB4\_Exp7\_Scramble en `ingest_stats.tsv` | 22.207.504 |
| A3 | Suma de `counts_mirtop.tsv` en NB4\_Exp7\_Scramble | No mayor que A2 y dentro de ±15 % de los 20.432.348 reads miRNA de mirTrace (es decir, entre 17,4 M y 22,2 M) |
| A4 | Secuencias multi-asignadas en el joined sin filtrar | 5.202 de 182.547 |
| A5 | Validación de entradas | 49 muestras, 10 contrastes, 4 `norm_group`, 0 avisos de validación |
| A6 | Flags de QC | LOW\_MIRNA\_FRAC = NB4\_Exp4\_DDX41sh1, NB4\_Exp6\_DDX41sh1, NB4\_Exp7\_DDX41sh1, MSCline\_Exp10\_DDX41sh1 y DDX41Patient\_MNC\_P6. HIGH\_TRNA\_FRAC = DDX41Patient\_MNC\_P6. HIGH\_RRNA\_FRAC = DDX41Patient\_MNC\_P3. LOW\_DEPTH: a lo sumo NB4\_Exp7\_DDX41sh1 (≈ 0,6 M reads miRNA según mirTrace; **\[VERIFICAR\]**) |
| A7 | Sentido del log2FC | En un conjunto sintético con un miRNA ×4 en treatment: log2FC\_mle ≈ +2 y `direction` = up |
| A8 | Reproducibilidad | Dos ejecuciones con la misma semilla dan tablas idénticas (md5) |
| A9 | DE | 10 contrastes × 3 ejecuciones producen `de_results.tsv` con padj no nulo; `de_summary.tsv` con 30 filas menos las omisiones documentadas |
| A10 | Fase 2: TRIM\_COLLAPSE frente a nf-core | Reads tras fastp dentro de ±3 % del Passed Filter de nf-core en las 49 muestras |
| A11 | Fase 2: M6 | Con la guía TTGATGACGTGCTGGATGTTT en `shrna_fasta`: fracción mayor que 0,85 en NB4\_Exp4\_DDX41sh1 y menor que 0,05 en cada Scramble |
| A12 | Fase 2: M7 | Suma de clases + unaligned + multimap\_gt\_max igual a los reads tras fastp (±0,5 %) en las 49 muestras |

### 9.2 Tests técnicos

- **nf-test** con un conjunto de prueba de 6 muestras (NB4\_Exp4, Exp6 y Exp7 de Scramble y de DDX41sh1; un solo contraste), generado por `bin/make_test_data.py` a partir del run real.
- **Pruebas unitarias** (`testthat` y `pytest`) del reparto 1/k de multi-asignadas, de los filtros de isomiR, del sentido del log2FC y del parser de metadata.
- **Casos de fallo** que deben terminar con un mensaje claro: muestra ausente en los conteos, contraste con un nivel inexistente, fichero obligatorio ausente, nombre de muestra que no casa con las reglas de `make_metadata.py`.

## 10. Convenciones y reglas para quien implemente

- Cada proceso lleva `tag` con el id de muestra o de grupo, `label` de recursos, `stub:` para pruebas y emite un `versions.yml`. Los `publishDir` se definen solo en `conf/modules.config`.
- Scripts en `bin/` con `--help`; R con `optparse`, Python con `argparse`; ninguna ruta absoluta dentro del código; `set.seed(params.seed)`; `sessionInfo()` volcado en `pipeline_info/`.
- Sin acceso a red en ejecución (todo va en el contenedor o en recursos locales). multiMiR (online) no se usa por defecto.
- Nunca se modifica la carpeta de nf-core. Nunca se descarta una muestra sin dejar rastro en `qc_flags.tsv` y en el informe.
- Todos los umbrales viven en `nextflow.config` y se vuelcan a `pipeline_info/params.json`. Un umbral no puede estar escrito dentro de un script.
- Mensajes de error con el patrón: qué falló, dónde, qué se esperaba y cómo arreglarlo.
- Las tablas de salida tienen un contrato fijo de columnas (secciones 7 y 8); cambiar columnas exige subir la versión del pipeline.
- No pasar a la fase N+1 hasta que los tests de aceptación de la fase N pasen.
- Cuando un punto marcado **\[VERIFICAR\]** no coincide con la realidad: parar, documentar lo encontrado y consultar.

## 11. Hitos de implementación

| Hito | Contenido | Definición de hecho |
| --- | --- | --- |
| 1 | Esqueleto, esquema de parámetros, perfiles y M0 | `nextflow run . -profile test,singularity` valida las entradas; A5 pasa |
| 2 | M1, M2 y M3 | A1-A4 y A6 pasan; figuras de QC y de composición generadas |
| 3 | M4, M5 y M10 (informe v1) | A7-A9 pasan; informe HTML con QC, PCA y DE |
| 4 | TRIM\_COLLAPSE y M6 | A10 y A11 pasan |
| 5 | M7 | A12 pasa; composición completa y Δ por contraste |
| 6 | M8 y M9 | Tablas y figuras de la fase 3; informe v2 |
| 7 | M11 (opcional) | Lista de candidatos novel con evidencia |

## 12. Cuestiones abiertas (necesitan un dato o una decisión)

1. `params.yaml` y `conf/custom.config` del run de nf-core (para alinear executor, recursos y parámetros reales).
2. Cabecera y contenido de `samplemap_rnaseq_2026_04.tsv` (fuente de `make_metadata.py`).
3. `generate_samplesheet.sh` apunta a `raw_data/snoRNAseq/...` pero el `samplesheet.csv` actual apunta a `raw_data/smRNAseq/...`: confirmar cuál es el correcto y regenerar si procede.
4. Secuencias de los shRNA (guía y, si es posible, la diana en el mRNA) de DDX41 sh1 y sh2, DHX34 sh1 y sh2 y Scramble. Hipótesis sobre sh1: guía TTGATGACGTGCTGGATGTTT (diana AAACATCCAGCACGTCATCAA).
5. Kit de librería usado por Novogene (el adaptador presente es el universal de Illumina).
6. Nombres reales de los ficheros de mirTrace y de MultiQC (`ls mirtrace/MUESTRA`, `ls multiqc/multiqc_data`).
7. Ubicación del índice Bowtie1 del genoma y release del GTF de Ensembl; fuente de las anotaciones de tRNA y de miRBase GFF3 (misma release que la usada en nf-core).
8. Decisión sobre relanzar nf-core con `--save_intermediates` (no hace falta para el MVP; solo si se quieren BAM o FASTQ recortados de nf-core).
9. Validar `consistency_pairs.csv` (sección 7, M5).
10. Interpretación biológica: si un shRNA domina la librería, el efecto sobre los miRNAs endógenos puede deberse en parte a una saturación de la maquinaria de RNAi y no solo al gen silenciado. Por eso la fracción shRNA y la concordancia sh1 frente a sh2 acompañan siempre a los resultados de DE.

## 13. Glosario

- **miRNA canónico:** secuencia exacta del miRNA maduro en miRBase. **isomiR:** variante con recortes, extensiones, adiciones no plantilladas o SNV respecto al canónico.
- **norm\_group:** subconjunto de muestras que comparten modelo de DE y normalización (aquí, un sistema celular o un tejido).
- **flag:** aviso de QC sobre una muestra; no implica exclusión.
- **multi-asignada:** secuencia que mirtop asigna a más de un miRNA; se reparte 1/k.
- **reference y treatment:** nivel base y nivel comparado de un contraste; log2FC = treatment frente a reference.
- **VST:** transformación estabilizadora de la varianza de DESeq2, usada para PCA y heatmaps (no para DE).
- **CLR:** transformación log-ratio centrada para datos composicionales (porcentajes por clase).
- **shRNA guía y pasajera:** las dos hebras del shRNA procesado; la guía es la complementaria a la diana.
