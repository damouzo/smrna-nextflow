# smrna_nextflow — Plan v3 (simplificado)

*Pipeline downstream de small RNA-seq sobre los resultados de nf-core/smrnaseq 2.4.x. Proyecto CRUK_DDX41_DHX34, batch_2026_04. Pocos scripts, cada uno con una función clara. Sin tests: se validará con ejecuciones reales.*

## 1. Alcance

**Contrato:** entrada = carpeta `results/` de nf-core/smrnaseq (`nfcore_dir`) + `metadata.csv` + `contrasts.csv`; salida = carpeta `results/` (`outdir`, por defecto dentro de la carpeta del pipeline `smrna_nextflow/`) con tres ramas: `deliverables/` (tablas y figuras para compartir con el wet lab), `raw/` (volcados nativos para auditoría interna) y `run_info/` (parámetros, umbrales y metadatos del run).

| Fase | Qué hace |
| --- | --- |
| 1 (MVP) | Validar entradas, matriz de conteos de miRNA, QC y flags, composición por biotipo (mirTrace), PCA, DE por contraste (DESeq2) con análisis de sensibilidad, informe |
| 2 | Composición completa por biotipo, con DE de snoRNA, tRNA, rRNA y snRNA |
| 3 | Perfil de isomiRs; enriquecimiento funcional con miRTarBase |

Fuera de alcance: repetir el alineamiento de miRNA (lo hace nf-core), miRNAs novel (miRDeep2), piRNA, predicción de dianas por secuencia.

Todas las muestras se tratan igual: no hay lógica específica para ningún tipo de experimento. Las diferencias entre muestras entran solo por `metadata.csv`, `contrasts.csv` y los flags genéricos de QC.

Reglas de diseño:

- Nada específico del batch en el código.
- Los datos de nf-core son de solo lectura. Todo se escribe en `outdir`.
- Sin CI/CD de tests, sin Seqera/Tower y sin plantilla nf-core. La única automatización es la GitHub Action que construye el contenedor.
- Hay internet en el clúster, también en los nodos de cómputo **[verificado con `srun`: mirbase.org y ghcr.io responden]**. Aun así, el pipeline no descarga nada durante el análisis: las referencias grandes se descargan una vez y se pasan por ruta (D8).

## 2. Hechos del batch_2026_04

Marcados **[verificado]** los comprobados sobre los datos al revisar este plan. El resto vienen del autor original.

| # | Hecho | Consecuencia |
| --- | --- | --- |
| H1 | 49 muestras single-end, 37-90 M reads tras fastp. Librería de ligación 3'/5' con inserto de 18-40 nt. Novogene da el adaptador 3' `AGATCGGAAGAGCACACGTCT` y el 5' `GTTCAGAGTTCTACAGTCCGACGATC` **[verificado]**. nf-core se lanzó con el perfil `illumina`, que declara `three_prime_adapter: TGGAATTCTCGGGTGCCAAGG`. En 100.000 reads crudos de C1, `AGATCGGAAGAGC` aparece en el 98,2 % y `TGGAATTCTCGGGTG` en el 0 % **[verificado]**. La página de servicio de Novogene no detalla kit ni adaptadores (solo SE 50 bp y 10-20 M reads recomendados). | Adaptador correcto: `AGATCGGAAGAGCACACGTCT`; el valor por defecto de nf-core (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCA`) recorta igual **[verificado]**. El FASTQ usado es `01.RawData` (con adaptador). MultiQC mostraba los adaptadores como eliminados porque FastQC busca el motivo universal de 12 nt y fastp sí recortaba con `known_adapters.fa` (ver H2); `fastp-pct_adapter` (0,002-0,07 %) solo cuenta lo recortado con el adaptador declarado. |
| H2 | **El recorte de nf-core dejó una `A` final en el 98,8 % de los reads** (200.000 reads de la muestra C1) **[verificado]**. Con el adaptador de Novogene (21 nt) o con el por defecto de nf-core (34 nt), en el mismo fastp y con o sin `known_adapters.fa`, el último nt es `A` en el 20,1 %, lo esperable por azar **[verificado]**. Causa: `known_adapters.fa` de fastp trae `GATCGGAAGAGCACACGTCTG`, sin la `A` inicial, y es lo único que recorta cuando el adaptador declarado no aparece. | En mirtop, el 62 % de los reads de NB4_Exp7_Scramble tiene `iso_add3p = 1` y solo el 2 % es canónico **[verificado]**. Explica el 9,6 % de `mature_counts` frente a mirtop. Los conteos de miRNA por mirtop siguen siendo válidos, pero `mature_counts` no es utilizable y el perfil de isomiRs 3' solo tiene sentido con lecturas bien recortadas. Ver D1. |
| H3 | 5.202 de 182.547 secuencias únicas (2,8 %) están asignadas a más de un miRNA. `mirna.tsv` suma esas filas, así que las cuenta dos veces. | Los conteos se calculan desde `joined_samples_mirtop.tsv`, repartiendo cada secuencia 1/k entre sus k miRNAs. Contando cada UID una vez, NB4_Exp7_Scramble suma 22,17 M **[verificado]**. |
| H4 | mirTrace: miRNA 1,5-79 % de los reads (media 44 %), tRNA hasta 52 %, rRNA hasta 33 %, unknown 19-98 %. | La composición por biotipo es muy variable entre muestras. |
| H5 | **Existen secuencias colapsadas y recortadas por muestra**: `mirtrace/<muestra>/qc_passed_reads.all.collapsed/<muestra>.fastp.fasta`, con cabecera `>seq_N_xCOUNT rnatype:CLASE` (5,4 GB en total) **[verificado]**. No hay BAM ni FASTQ recortados. | M6 lee estos FASTA. No hace falta FASTQ crudo ni re-recortar. |
| H6 | `mature_counts.csv` está traspuesto (muestras en filas, 2.301 miRNAs hsa) **[verificado]**. `mirna.tsv` y `joined_samples_mirtop.tsv` traen las muestras en orden arbitrario. | Unir siempre por nombre de muestra. |
| H7 | `mature` y `hairpin` del run apuntan a `nf-core/test-datasets`: es miRBase 22.1 (2.656 miRNAs hsa, idéntico a `mirbase.org/download_version_files/22.1/mature.fa`). Los 2.301 de `mature_counts.csv` son los detectados **[verificado]**. `mirna_gtf` era null y nf-core descargó el GFF3 `CURRENT` de miRBase de ese día (la cabecera `miRBasev22` de los GFF de mirtop es un valor por defecto de la herramienta y no prueba la versión). Hoy el miRBase vigente es la 23.1 (sep 2026), con 20.514 nombres aclarados y 9.091 secuencias maduras ajustadas **[verificado]**. | Un relanzamiento con `mirna_gtf` vacío mezclaría mature y hairpin de la 22.1 con un GFF3 de la 23.1. Fijar los tres a la misma versión (D10). |
| H8 | Hay 2 sistemas × 5 condiciones × 3 experimentos (NB4: Exp4, 6, 7; MSCline: Exp8, 9, 10), cada experimento con todas las condiciones **[verificado]**. | El diseño es pareado por experimento. Ver D5. |
| H9 | El proyecto se pidió a Novogene como "snoRNAseq", pero la librería es de small RNA con inserto de 18-40 nt **[verificado]**. | Los snoRNA completos (60-300 nt) no se capturan: solo fragmentos derivados. Se avisa en el informe. |
| H10 | mirTrace descarta los reads de menos de 18 nt (`LENGTH_SHORTER_THAN_18`: 1,5 M reads en NB4_Exp7_Scramble) **[verificado]**, mientras fastp conservó los de 17 nt (`fastp_min_length: 17`). | Los FASTA colapsados (H5) no incluyen reads de 17 nt. Con `fastp_min_length: 18` los totales de nf-core y de M6 coinciden. |

## 3. Decisiones de diseño

| ID | Decisión | Motivo |
| --- | --- | --- |
| D1 | **Relanzar nf-core 2.4.1** (decidido) con los cambios de la sección 4.5. Sin relanzar, el DE seguiría siendo válido, pero se desactivaría el perfil de isomiRs 3' y M6 usaría `trim3p = 1`. | H1, H2, H7, H10 |
| D2 | Nextflow DSL2 plano: `main.nf`, `nextflow.config`, `modules/`. Un script de R por etapa en `bin/`. Parámetros en `nextflow.config`, sin esquema nf-schema. | Simplicidad |
| D3 | **Un único contenedor** `ghcr.io/<owner>/smrna-nextflow` (R/Bioconductor + bowtie + bedtools), construido por GitHub Actions desde `containers/smrna/Dockerfile`, con el mismo patrón que `lncRNA-nextflow/.github/workflows/build-containers.yml`. Sin conda, sin renv. | Un solo sitio donde cambian las versiones |
| D4 | `count_source = mirtop` por defecto; `mature` solo como sensibilidad. Reparto 1/k de multi-asignadas y redondeo final a entero (fijos). | H2, H3 |
| D5 | Un objeto DESeq2 por `norm_group`, modelo principal `~ condition`; log2FC = treatment/reference. Si se indica `batch_col` (aquí `exp`), una ejecución extra de sensibilidad usa `~ batch_col + condition`. | H8; `exp` no entra en el modelo principal, decisión previa |
| D6 | Las muestras marcadas por QC no se excluyen; el DE se repite sin ellas como sensibilidad. | Evitar exclusiones silenciosas |
| D7 | Composición por biotipo en dos niveles: mirTrace (fase 1) y completa por alineamiento (fase 2). | Empezar simple |
| D8 | **Referencias grandes por ruta, nunca descargadas por el pipeline.** El README explica cómo descargar cada una con `wget`; el usuario pone la ruta en `params.yaml` o en `nextflow.config`. | Reproducibilidad y simplicidad |
| D9 | Semilla fija (`seed = 42`). Versiones y parámetros se vuelcan a `run_info/pipeline_info/`. | Reproducibilidad |
| D10 | **miRBase 23.1 fijada por ruta**, nunca por URL `CURRENT`: `mature`, `hairpin` y `mirna_gtf` apuntan a `data/reference/genome/GRCh38/mirbase/v23.1/` (ya descargada, con README y md5). Para actualizar en el futuro se descarga la nueva release en `mirbase/vX.Y/` y se cambian las tres rutas. La 22.1 no es viable como juego completo: no hay `hsa.gff3` de esa versión. | H7. Los nombres cambian respecto al run de abril y puede que respecto a miRTarBase (no verificado): M8 lista los miRNAs no emparejados. |

## 4. Entradas

### 4.1 Ficheros de nf-core que se leen

Rutas relativas a `nfcore_dir` (la carpeta `results/`).

| Fichero | Contenido | Lo usa |
| --- | --- | --- |
| `mirna_quant/mirtop/joined_samples_mirtop.tsv` | Columnas `UID, Read, miRNA, Variant, iso_5p, iso_3p, iso_add3p, iso_snp, …_nt` y una columna por muestra. Una fila por par (secuencia, miRNA). | M1, M7 |
| `mirna_quant/mirtop/mirna.tsv` | miRNA × muestra, con doble conteo | solo cruce en M1 |
| `mirna_quant/edger_qc/mature_counts.csv` | muestras × miRNA (traspuesto) | M1 (sensibilidad) |
| `multiqc/multiqc_data/multiqc_general_stats.txt` | estadísticas de fastp y mirtop **[verificado]** | M1 |
| `multiqc/multiqc_data/multiqc_samtools_stats.txt` | % mapeado al genoma | M1 |
| `mirtrace/<muestra>/mirtrace-stats-rnatype.tsv` | categorías en filas (`RNA_TYPE`) y una columna de muestra **[verificado]** | M1, M2 |
| `mirtrace/<muestra>/qc_passed_reads.all.collapsed/<muestra>.fastp.fasta` | secuencias únicas con conteo y clase (H5) | M0, M6 |
| `pipeline_info/params_*.json` | parámetros del run | M0 |

No se usan `mirdeep2/`, `genome_quant/bam/` (solo trae `.stats`), `mirna_quant/mirtop/export/` ni los FASTQ.

### 4.2 `metadata.csv`

Una fila por muestra. Las muestras sin librería se omiten (aquí DDX41Patient_MSC_P6).

| Columna | Obligatoria | Descripción |
| --- | --- | --- |
| `sample_id` | sí | Igual al nombre de muestra de nf-core |
| `condition` | sí | Etiqueta usada en `contrasts.csv` |
| `norm_group` | sí | Subconjunto con modelo y normalización propios (NB4, MSCline, MNC, MSC) |
| `exp` | no | Experimento (Exp4…); color de PCA y `batch_col` |
| `subject` | no | Donante o paciente |
| `library_id` | no | C1…C49 |
| `wts_id` | no | Identificador de RNA-seq |

`assets/metadata.csv` es una plantilla genérica. La del batch está en `demo/batch_2026_04/metadata.csv`, generada con `demo/batch_2026_04/make_metadata.R` a partir de `samplemap_rnaseq_2026_04.tsv`. Ese script es de un solo uso y no forma parte del pipeline.

### 4.3 `contrasts.csv`

```
contrast_id,variable,reference,treatment,norm_group
NB4_DDX41_sh1_vs_Scramble,condition,NB4_Scramble,NB4_DDX41_sh1,NB4
MNC_DDX41_Patient_vs_Healthy,condition,MNC_Healthy,MNC_DDX41_Patient,MNC
```

`demo/batch_2026_04/contrasts.csv` lleva los 10 contrastes del batch. M0 comprueba que `reference` y `treatment` existen en `variable`, pertenecen al `norm_group`, difieren y tienen al menos 2 muestras por nivel (aviso con menos de 3). log2FC positivo significa más abundante en `treatment`.

### 4.4 Referencias (D8)

Solo las necesita la fase 2 y la fase 3. Se pasan por ruta. En este proyecto ya están casi todas en `data/reference/genome/GRCh38/` (`demo/batch_2026_04/params.yaml` las referencia); el README trae los `wget` para quien parta de cero.

| Recurso | Parámetro | Estado en este proyecto |
| --- | --- | --- |
| Genoma Ensembl (`dna_sm.primary_assembly.fa.gz`) | `genome_fasta` | existe; contigs `1`, `MT` (sin `chr`) |
| GTF Ensembl 115 | `gtf` | existe; solo trae 22 `Mt_tRNA` (sin tRNA nucleares) **[verificado]** |
| Anotación de ncRNA de RNAcentral (`homo_sapiens.GRCh38.bed`) | `ncrna_bed` | **ya descargado: se reutiliza el mismo fichero de lncRNA-nextflow**, `data/reference/genome/GRCh38/rnacentral/homo_sapiens.GRCh38.bed` (BED12, 21-oct-2026). Tiene 1.016 contigs `chr…` con alt y random; se usan solo `chr1-22,X,Y,M`, porque el FASTA primary assembly no trae los alt. Clase en la columna 14: tRNA 7.370, snoRNA 2.679, rRNA 53.707, miRNA 6.978, snRNA 3.110, Y_RNA 1.038 **[verificado]**. Trae también circRNA, piRNA y lncRNA, que se descartan al construir la anotación. |
| Índice Bowtie1 del genoma | `bowtie_index` | no existe; si es null, M6 lo construye una vez desde `genome_fasta` |
| miRTarBase (MTI de la especie) | `mirtarbase_file` | descargar a mano (fase 3) |

### 4.5 Relanzamiento de nf-core: qué cambiar

miRBase 23.1 ya está descargada y fijada (D10). Las URL `CURRENT` cambian con cada release, por eso la versión va en el nombre de la carpeta. `mirbase.org` sirve `hsa.gff3` como una sola línea envuelta en HTML, así que se convierte a GFF3 plano (el original queda como `hsa.gff3.html_orig`); los tres ficheros son coherentes entre sí (mismos hairpins, mismos miRNAs maduros salvo 6 sin coordenadas) **[verificado]**:

```bash
sed 's/<br>/\n/g; s/<\/\?p>//g; s/&quot;/"/g; s/&amp;/\&/g' hsa.gff3.html_orig > hsa.gff3
```

Los ficheros finales están en `analysis/smRNAseq/batch_2026_04/nfcore_smrnaseq/`: `params.yaml` y `run_command.sh` (el anterior queda como `run_command_v1_april.sh`). `-preview` de Nextflow valida los parámetros y el samplesheet **[verificado]**. Contenido relevante de `params.yaml`:

```yaml
three_prime_adapter: "AGATCGGAAGAGCACACGTCT"   # adaptador de Novogene (H1)
fastp_min_length: 18                          # coincide con el mínimo de mirTrace (H10)
mature:    /data/BCI-KRP/projects/CRUK_DDX41_DHX34/data/reference/genome/GRCh38/mirbase/v23.1/mature.fa
hairpin:   /data/BCI-KRP/projects/CRUK_DDX41_DHX34/data/reference/genome/GRCh38/mirbase/v23.1/hairpin.fa
mirna_gtf: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/data/reference/genome/GRCh38/mirbase/v23.1/hsa.gff3
skip_mirdeep: true                            # opcional: el pipeline no usa miRDeep2
```

`run_command.sh`: `-profile singularity` sin `illumina`, `--outdir results` (el run de abril queda intacto en `results/`) y `-work-dir /gpfs/scratch/qp241615/smrna-cruk-batch_2026_04`.

No cambiar: `filter_contamination: false` (con true, nf-core retiraría rRNA, tRNA y otros ncRNA antes de cuantificar y se perdería la composición), `fastp_max_length: 100`, `mirtrace_species: hsa`, `with_umi: false`. `mature.fa` trae todas las especies y nf-core filtra `hsa` con `seqkit` **[verificado]**.

Qué aporta cada pieza: nf-core cuantifica solo miRNA (mirtop). snoRNA, tRNA, rRNA y el resto de especies salen de los FASTA colapsados de mirTrace (H5), que nf-core publica por defecto, mediante M6. `save_intermediates: true` solo añade los BAM de genoma, que este pipeline no usa.

## 5. Parámetros

Todos tienen valor por defecto en `nextflow.config` y se vuelcan a `run_info/pipeline_info/params.json`. Un umbral no puede estar escrito dentro de un script.

**Por ejecución** (`params.yaml`):

| Parámetro | Default | Descripción |
| --- | --- | --- |
| `nfcore_dir`, `metadata`, `contrasts` | obligatorios | Entradas |
| `outdir` | `results` | Salida; contiene `deliverables/`, `raw/` y `run_info/` |
| `batch_col` | null | Columna de `metadata` para la ejecución `sens_batch` |
| `run_biotype`, `run_isomir`, `run_functional` | false | Activan las fases 2 y 3 |

**Por proyecto**:

| Parámetro | Default | Descripción |
| --- | --- | --- |
| `species` | `hsa` | Prefijo de miRNA |
| `expected_adapter` | `AGATCGGAAGAGCACACGTCT` | M0 avisa si el run declaró otro |
| `trim3p` | 0 | Nt a recortar del extremo 3' en M6 (1 si no se relanza nf-core, D1) |
| `genome_fasta`, `bowtie_index`, `gtf`, `ncrna_bed` | null | Referencias de M6 |
| `mirtarbase_file` | null | M8 |
| `container` | null | Sustituye la imagen por defecto (D3) |

**Análisis** (defaults, se tocan solo por decisión de análisis):

| Parámetro | Default | Descripción |
| --- | --- | --- |
| `count_source` | `mirtop` | `mirtop` o `mature` |
| `isomir_exclude_variants` | `iso_snv_seed,iso_snv_central,iso_snv_central_offset,iso_snv_central_supp` | Filas con estas etiquetas en `Variant` se descartan |
| `isomir_max_shift_5p`, `isomir_max_shift_3p`, `isomir_min_len` | 2, 5, 18 | Filtros de isomiR |
| `qc_min_mirna_frac`, `qc_max_trna_frac`, `qc_max_rrna_frac`, `qc_min_mirna_reads` | 0.10, 0.30, 0.20, 1000000 | Umbrales de flags |
| `sensitivity_exclude_flags` | `LOW_MIRNA_FRAC,LOW_DEPTH` | Flags que excluye `sens_no_flagged` |
| `min_cpm`, `min_total_count` | 1, 10 | Filtro de expresión; `min_samples` = tamaño del nivel más pequeño del grupo |
| `alpha`, `lfc_flag` | 0.05, 1 | FDR y umbral de log2FC para "strong" |
| `shrink_type`, `n_top_var` | `ashr`, 500 | Contracción de log2FC y miRNAs para la PCA |
| `bowtie_mismatches`, `bowtie_max_hits`, `min_overlap_frac` | 1, 50, 0.8 | M6 |
| `biotype_priority` | `miRNA,tRNA,rRNA,snoRNA,snRNA,scaRNA,Y_RNA,vault_RNA,lncRNA,protein_coding,unannotated` | Orden de la asignación jerárquica |
| `de_classes` | `snoRNA,snRNA,tRNA,rRNA` | Clases con DE a nivel de feature |
| `isomir_min_reads`, `min_sig_mirnas` | 500, 5 | M7 y M8 |
| `seed` | 42 | Semilla |

Las figuras se guardan siempre en PDF y PNG.

## 6. Estructura del repositorio

```
smrna-nextflow/
├── README.md               # corto y limpio (ver sección 11)
├── main.nf                 # el flujo completo, legible de arriba abajo
├── nextflow.config         # parámetros, imagen, perfiles (singularity, docker, slurm)
├── modules/
│   ├── inputs.nf           # M0 y M1
│   ├── analysis.nf         # M2, M3 y M4
│   ├── reads.nf            # M6 (fase 2)
│   └── extras.nf           # M7 y M8 (fase 3), M5 informe
├── bin/                    # un script de R por etapa, con --help (optparse)
│   ├── validate_inputs.R
│   ├── build_counts.R
│   ├── qc_composition.R
│   ├── deseq2_group.R
│   ├── summarize_de.R
│   ├── biotype_full.R      # fase 2
│   ├── isomir_profile.R    # fase 3
│   ├── functional.R        # fase 3
│   └── report.Rmd
├── containers/smrna/Dockerfile
├── .github/workflows/build-containers.yml
├── assets/                 # plantillas genéricas: metadata.csv, contrasts.csv, params.yaml
├── demo/batch_2026_04/     # ejemplo real: metadata, contrasts, params.yaml, make_metadata.R, run.sh
└── docs/
```

Contenedor (D3): la imagen parte de `bioconductor/bioconductor_docker` con DESeq2, ashr, ggplot2, ggrepel, pheatmap, data.table, optparse, jsonlite, yaml, rmarkdown y clusterProfiler, más `bowtie` y `bedtools` de apt. `nextflow.config` asigna una sola imagen por defecto a todos los procesos, con `params.container` para sustituirla. GitHub Actions la construye y publica en GHCR al hacer push a `main` si cambia algo en `containers/`, o a mano con `workflow_dispatch`. Para que Singularity la descargue en el clúster, el paquete de GHCR debe ser público.

Perfiles: `singularity` y `docker` (motor) y `slurm` (executor, cola `compute`, máximos de 32 CPU, 128 GB y 72 h, como `conf/custom.config` del run de nf-core). `demo/batch_2026_04/run.sh` replica el patrón de `run_command.sh` de nf-core: módulo de nextflow, `NXF_VER=25.04.3`, `-work-dir` en `/gpfs/scratch/…` (nunca dentro del proyecto) y `-resume`.

Reglas: ninguna ruta absoluta dentro del código; `set.seed(params.seed)`; cada proceso con `tag`, `label` de recursos y `publishDir` definido junto al proceso; mensajes de error con qué falló, dónde, qué se esperaba y cómo arreglarlo.

## 7. Flujo

```
Fase 1:  M0 validar ─> M1 conteos ─> M2 QC y composición ─> M3 PCA y DE ─> M4 resumen ─> M5 informe
Fase 2:                              M6 biotipo completo   (lee los FASTA colapsados de mirtrace)
Fase 3:                              M7 isomiRs   M8 funcional
```

## 8. Módulos

### M0 · `validate_inputs.R` (fase 1)

- **Entrada:** `nfcore_dir`, `metadata`, `contrasts`.
- **Hace:**
  1. Comprobar que existen los ficheros de la tabla 4.1.
  2. Leer `pipeline_info/params_*.json` y avisar (no fallar) si el adaptador declarado difiere de `expected_adapter`, si `with_umi` o `filter_contamination` son true, o si `mature`, `hairpin` y `mirna_gtf` no están en la misma carpeta (o son URL o `test-datasets`). Las tres rutas, que llevan la versión de miRBase, se anotan en `audit.md`.
  3. Mirar el último nt de los reads de un FASTA colapsado; si más del 90 % termina en la misma base, avisar de resto de adaptador (H2).
  4. Validar `metadata.csv` (columnas obligatorias, `sample_id` únicos, todas las muestras presentes en los conteos y ninguna sobrante) y `contrasts.csv` (sección 4.3).
- **Salida** (`run_info/`): `audit.md` (versión de nf-core, genoma, adaptador, avisos, nº de muestras por condición y grupo), `metadata_validated.csv`, `contrasts_validated.csv`.
- **Falla** si faltan ficheros, no casan las muestras o un contraste es inválido.

### M1 · `build_counts.R` (fase 1)

- **Entrada:** ficheros de 4.1 y `metadata_validated.csv`.
- **Hace (fuente mirtop):**
  1. Leer `joined_samples_mirtop.tsv` con `data.table::fread` y quedarse con las muestras de `metadata`.
  2. Descartar filas con alguna etiqueta de `isomir_exclude_variants` en `Variant`, con `|iso_5p|` o `|iso_3p|` por encima de los máximos, o con `Read` más corto que `isomir_min_len`. Contar los reads retirados por muestra.
  3. k = nº de miRNAs distintos por UID tras el filtrado; peso 1/k; conteo(muestra, miRNA) = suma de peso × conteo; redondeo a entero.
- **Hace (fuente mature):** leer `mature_counts.csv`, trasponer y ordenar por `metadata`.
- **Qc table:** fusiona fastp, mirTrace, samtools, totales de mirtop y metadata.
- **Salida:** `deliverables/counts/counts_mirtop.tsv` (miRNA × muestra, enteros, columnas en el orden de `metadata`) y `deliverables/qc/qc_table.tsv` con `sample_id, condition, norm_group, exp, reads_raw, reads_passed_fastp, pct_passed, mirtrace_pct_miRNA, mirtrace_pct_tRNA, mirtrace_pct_rRNA, mirtrace_pct_artifact, mirtrace_pct_unknown, mirtop_raw_total, mirna_reads_mirtop, mature_total, canonical_frac, multi_assigned_frac, genome_mapped_pct`. Los volcados de sensibilidad van a `raw/counts/`: `counts_mature.tsv` e `isomir_filtered.tsv.gz` (filas supervivientes con su peso).
- Los miRNAs de `mature` que no casan por nombre con mirtop se listan en el log, no son un error.

### M2 · `qc_composition.R` (fase 1)

- **Entrada:** `qc_table.tsv`, metadata y contrasts.
- **Flags** (se marcan, no se excluyen): LOW_MIRNA_FRAC, HIGH_TRNA_FRAC, HIGH_RRNA_FRAC y LOW_DEPTH.
- **Composición mirTrace:** % de miRNA, rRNA, tRNA, artifact y unknown por muestra; por contraste, media por condición y Δ en puntos porcentuales (descriptivo, sin test).
- **Salida** (`deliverables/qc/`): `qc_flags.tsv` (`sample_id, flag, value, threshold`), `flagged_samples.txt`, `qc_overview.pdf/png`, `composition_by_sample.tsv`, `composition_delta.tsv` (`contrast_id, category, mean_ref, mean_trt, delta_pp`), `composition_stacked_bar.pdf`.

### M3 · `deseq2_group.R` (fase 1)

Un proceso por `norm_group` (más uno `all` sin DE con todas las muestras) y ejecución.

- **Ejecuciones:** `main` (fuente `count_source`, todas las muestras), `sens_mature`, `sens_no_flagged` y, si hay `batch_col`, `sens_batch` (`~ batch_col + condition`). Un contraste con menos de 2 muestras en un nivel se omite y se registra.
- **Hace:**
  1. Filtro de expresión (CPM ≥ `min_cpm` en ≥ `min_samples`, suma ≥ `min_total_count`).
  2. DESeqDataSet; size factors con `ratio` y paso automático a `poscounts` si falla; `DESeq` con `fitType` paramétrico y paso a `local` o `mean` si falla (se registra).
  3. Solo en `main`: VST con `blind = TRUE`, PCA con los `n_top_var` miRNAs más variables (color `condition`, forma `exp`, muestras con flag resaltadas) y mapa de calor de distancias entre muestras.
  4. Por contraste: `results()` y `lfcShrink(type = shrink_type)`; `direction` = up, down o ns (por `alpha` y signo del log2FC); `strong` si |log2FC| ≥ `lfc_flag`.
- **Salida:** `deliverables/explore/GRUPO/` (`pca.pdf`, `pca_coords.tsv`, `sample_distance_heatmap.pdf`); `deliverables/de/CONTRAST_ID/` (`de_full.tsv`, `de_significant.tsv`, `volcano.pdf`) para `main`. Las sensibilidades van a `raw/deseq2/RUN/CONTRAST_ID/de_full.tsv` y `de_significant.tsv`; `size_factors.tsv`, `filtered_counts.tsv` y `dds.rds` a `raw/deseq2/GRUPO/`.
- **Columnas de `de_full.tsv`:** `mirna, baseMean, log2FC_mle, log2FC_shrunk, lfcSE, stat, pvalue, padj, direction, strong, n_ref, n_trt, mean_norm_count_ref, mean_norm_count_trt`. `de_significant.tsv` es el subconjunto con `padj < alpha` (mismas columnas).

### M4 · `summarize_de.R` (fase 1)

- **Salida** (`deliverables/de/`): `de_summary.tsv` (`run, contrast_id, n_tested, n_up, n_down, n_strong`) y `concordance.tsv` (`contrast_id, run_a, run_b, n_common, spearman_log2fc, n_sig_shared`), `main` frente a cada sensibilidad.

### M5 · `report.Rmd` (fase 1, crece por fases)

Informe HTML que solo lee `outdir` (`deliverables/` y `run_info/`) y falla si falta un fichero de las fases activas. Secciones: resumen, auditoría, QC y flags, composición, PCA, DE por contraste con concordancia, isomiRs y funcional (si existen), métodos y limitaciones.

Limitaciones (texto fijo): los porcentajes de composición son relativos; el DE es relativo al total de miRNA, no absoluto; solo se capturan fragmentos de snoRNA (H9); n pequeño.

**Salida:** `deliverables/report/smrna_report.html`.

### M6 · `biotype_full.R` (fase 2, `run_biotype`)

- **Entrada:** FASTA colapsados de mirtrace (H5), `bowtie_index` o `genome_fasta`, `gtf`, `ncrna_bed`, metadata y contrasts. Los FASTA solo contienen reads de 18 nt o más (H10).
- **Hace:**
  1. Una vez, construir `annotation.bed` con las clases de `biotype_priority`: ncRNA desde `ncrna_bed` (descartando circRNA y piRNA) y `protein_coding` y `lncRNA` desde el GTF. Nombres de cromosoma armonizados (quitar `chr`, `chrM` a `MT`); los contigs sin equivalente en el FASTA se cuentan en el log.
  2. `bowtie -f -v bowtie_mismatches -k (bowtie_max_hits + 1) --best --strata -3 trim3p`. Sin alineamiento: `unaligned`; más de `bowtie_max_hits`: `multimap_gt_max`.
  3. `bedtools intersect -wo -s -f min_overlap_frac` (misma hebra). Sin feature solapada: `unannotated`.
  4. Clase de mayor prioridad entre las solapadas; empate: `feature_id` alfabético. Una secuencia de tamaño n con N loci aporta n/N a cada locus.
  5. Para cada clase de `de_classes`, matriz feature × muestra (tRNA por anticodón) y DE con `deseq2_group.R`.
- **Salida** (`deliverables/biotype/`): `biotype_counts.tsv`, `biotype_composition.tsv`, `composition_delta.tsv`, `composition_stacked_bar.pdf` y `de/CLASE/CONTRAST_ID/{de_full.tsv, de_significant.tsv}`. Los volcados nativos `alignment_stats.tsv` y `feature_counts/CLASE.tsv` van a `raw/biotype/`.
- **Advertencia para el informe:** los porcentajes son composicionales; interpretar junto a los conteos absolutos. Los DE de clases con pocas features son exploratorios.

### M7 · `isomir_profile.R` (fase 3, `run_isomir`)

Un isomiR es un read de un miRNA que difiere de la secuencia canónica de miRBase: desplazamiento en 5' (`iso_5p`, cambia la semilla y por tanto las dianas), recorte o extensión en 3' (`iso_3p`), adición no plantillada en 3' (`iso_add3p`, normalmente A o U) o cambio interno (`iso_snp`).

- **Entrada:** `isomir_filtered.tsv.gz`, metadata, contrasts. **Requiere lecturas con el 3' bien recortado (D1).**
- **Hace:** por muestra y miRNA con al menos `isomir_min_reads`, fracciones canónica, con desplazamiento 5', con desplazamiento 3' y con adición 3', e isoforma dominante. Por contraste, cambio de `frac_5p_shifted` (t de Welch sobre logit, indicativo). Matriz de isomiRs (miRNA + iso_5p:iso_3p:iso_add3p) y DE con `deseq2_group.R`.
- **Salida** (`deliverables/isomir/`): `isomir_profile_by_sample.tsv`, `isomir_counts.tsv`, `de_isomir/CONTRAST_ID/{de_full.tsv, de_significant.tsv}`, `delta_5p_shift_CONTRAST_ID.tsv` y figuras.

### M8 · `functional.R` (fase 3, `run_functional`)

- **Entrada:** `de_full.tsv` de `main` y `mirtarbase_file`.
- **Hace:** targets validados de `species`, casados por nombre de miRNA. Por contraste y dirección, con al menos `min_sig_mirnas` significativos, enriquecimiento GO BP con clusterProfiler; el universo son los genes con algún target validado de un miRNA expresado.
- **Salida** (`deliverables/functional/CONTRAST_ID/`): `targets.tsv`, `ora_up.tsv`, `ora_down.tsv`.

## 9. Estructura de salida

Todos los productos viven bajo `outdir` (`results/` por defecto), organizados en tres ramas: `deliverables/` (lo que se comparte con el wet lab), `raw/` (volcados nativos para auditoría interna) y `run_info/` (parámetros y metadatos del run). En este proyecto el `outdir` vive en `analysis/smRNAseq/batch_2026_04/smrna_nextflow/`, aparte del run de nf-core (`nfcore_smrnaseq/results/`).

```
results/
├── deliverables/                        # tablas y figuras para el wet lab
│   ├── README.md                        # qué mirar y diccionario de columnas
│   ├── counts/      counts_mirtop.tsv
│   ├── qc/          qc_table.tsv, qc_flags.tsv, flagged_samples.txt,
│   │                composition_by_sample.tsv, composition_delta.tsv, plots/
│   ├── explore/     all/ y GRUPO/: pca.pdf, pca_coords.tsv, sample_distance_heatmap.pdf
│   ├── de/          de_summary.tsv, concordance.tsv
│   │                CONTRAST_ID/{de_full.tsv, de_significant.tsv, plots/volcano.pdf}
│   ├── biotype/     biotype_*.tsv, composition_*, de/CLASE/CONTRAST_ID/…              [fase 2]
│   ├── isomir/      isomir_profile_by_sample.tsv, isomir_counts.tsv,
│   │                de_isomir/CONTRAST_ID/…, delta_5p_shift_*.tsv, figuras           [fase 3]
│   ├── functional/  CONTRAST_ID/{targets.tsv, ora_up.tsv, ora_down.tsv}              [fase 3]
│   └── report/      smrna_report.html
├── raw/                                 # auditoría interna; no hace falta para interpretar
│   ├── README.md
│   ├── counts/      counts_mature.tsv, isomir_filtered.tsv.gz
│   ├── deseq2/      GRUPO/ (dds.rds, size_factors.tsv, filtered_counts.tsv)
│   │                sens_mature/ sens_no_flagged/ [sens_batch/]: CONTRAST_ID/{de_full.tsv, de_significant.tsv}
│   └── biotype/     alignment_stats.tsv, feature_counts/CLASE.tsv                     [fase 2]
└── run_info/                            # parámetros y metadatos del run
    ├── README.md
    ├── audit.md
    ├── metadata_validated.csv, contrasts_validated.csv
    ├── params_used.tsv, thresholds_used.tsv, run_manifest.tsv
    ├── sample_index.tsv, software_versions.tsv, qa_report.tsv
    ├── design_check/GRUPO.tsv
    └── pipeline_info/ params.json, software_versions.yml, execution_report.html, execution_trace.txt
```

| Pregunta | Archivo |
| --- | --- |
| ¿Están bien las muestras? | `deliverables/qc/qc_overview.pdf`, `deliverables/qc/qc_flags.tsv` |
| ¿Cómo es la composición del small RNA? | `deliverables/qc/composition_stacked_bar.pdf` (y `deliverables/biotype/…`) |
| ¿Se agrupan las muestras por biología? | `deliverables/explore/GRUPO/pca.pdf` |
| ¿Qué miRNAs cambian en cada contraste? | `deliverables/de/CONTRAST_ID/de_full.tsv` y `plots/volcano.pdf` |
| ¿Cuáles son significativos? | `deliverables/de/CONTRAST_ID/de_significant.tsv` |
| ¿Es robusto el resultado? | `deliverables/de/concordance.tsv` |
| Trazabilidad de parámetros y umbrales | `run_info/params_used.tsv`, `run_info/thresholds_used.tsv` |
| Resumen para compartir | `deliverables/report/smrna_report.html` |

Convenciones: snake_case; TSV con cabecera, tabuladores, sin comillas y con punto decimal; `sample_id` o `contrast_id` como primera columna cuando aplique; miRNAs con nombre completo (hsa-miR-21-5p). Cada contraste lleva `de_full.tsv` (todo lo probado) y `de_significant.tsv` (subconjunto con `padj < alpha`). Cambiar las columnas de una tabla de salida exige subir la versión del pipeline.

## 10. Hitos

| Hito | Contenido | Hecho cuando |
| --- | --- | --- |
| 1 | `main.nf`, `nextflow.config`, contenedor y Action, M0, M1 y M2 | El run sobre el batch valida las entradas y genera conteos, `qc_table.tsv`, flags y composición |
| 2 | M3, M4 y M5 | DE de los 10 contrastes en `main` y en las sensibilidades, e informe HTML con QC, PCA y DE |
| 3 | M6 | Composición completa y DE de las clases de `de_classes` |
| 4 | M7 y M8 | Tablas y figuras de la fase 3 |

No se pasa a la fase siguiente sin haber ejecutado la anterior sobre el batch real.

## 11. README

Un README corto y sin relleno: qué hace (2 líneas), requisitos (Nextflow, Singularity o Docker), inicio rápido (3 comandos), tabla de entradas, parámetros principales, **recursos opcionales con su `wget`** (D8), salidas y fases. Nada de historia del proyecto ni de decisiones; eso queda en `docs/`.

## 12. Cuestiones abiertas: resueltas y pendientes

| # | Cuestión | Resolución |
| --- | --- | --- |
| 1 | `params.yaml` y `custom.config` del run de nf-core | **Resuelta.** Executor slurm con `queueSize 50`; máximos de 32 CPU, 128 GB y 72 h; caché de Singularity en `…/CRUK_DDX41_DHX34/.singularity_cache`; `igenomes_ignore: true`; FASTA propio; `mirtrace_species: hsa`; `fastp_min_length 17`, `fastp_max_length 100`; `three_prime_adapter` del perfil `illumina` (incorrecto, H1); `with_umi` y `filter_contamination` en false; `mature`/`hairpin` de test-datasets (H7). |
| 2 | Cabecera del samplemap | **Resuelta.** `data/metadata/samplemap_rnaseq_2026_04.tsv`, 50 filas: `sample_id, model_system, cell_type, experiment, condition, donor_id, wts_rnaseq_sample_id, sno_rnaseq_sample_id`. DDX41Patient_MSC_P6 no tiene librería (NA), de ahí las 49 muestras. En el samplemap el sistema se escribe `MSC_line` y en los nombres de muestra `MSCline`. |
| 3 | Ruta del FASTQ (`snoRNAseq` frente a `smRNAseq`) | **Resuelta.** Existe `raw_data/smRNAseq/batch_2026_04/fastq/C1…C49`; `snoRNAseq` no existe. El `samplesheet.csv` actual ya apunta a `smRNAseq`; `generate_samplesheet.sh` está desfasado y conviene corregirlo antes de reutilizarlo. |
| 4 | Kit de librería | **Resuelta.** Ligación de adaptadores 3' y 5', RT y PCR, inserto de 18-40 nt (informe de Novogene). Adaptadores en H1. |
| 5 | Nombres de los ficheros de mirTrace y MultiQC | **Resuelta.** Ver tabla 4.1. |
| 6 | Índice Bowtie1, GTF, anotación de tRNA y miRBase | **Resuelta.** GTF Ensembl 115 y `rnacentral/homo_sapiens.GRCh38.bed` existen. El índice Bowtie1 no existe y lo construye M6. No hace falta miRBase GFF3 para este pipeline. |
| 7 | Relanzar nf-core | **Decidido** (D1, sección 4.5). El adaptador correcto elimina la `A` residual en una prueba con 200.000 reads (98,8 % → 20,1 %). |
| 8 | Internet en nodos de cómputo | **Resuelta.** Un `srun` en la cola `compute` alcanza mirbase.org y ghcr.io. |
| 9 | miRTarBase | **Pendiente de ti** (fase 3): descargarlo y pasar `mirtarbase_file`. |
| 10 | Versión de miRBase | **Resuelta.** 23.1 fijada por ruta (D10). |
