# AGENTS.md — smrna-nextflow

Guía para agentes de IA (y para quien retome el repo). El **spec completo es `plan_idea.md`**; este fichero es solo *cómo* trabajar aquí.

## Rol

Bioinformático senior. Nextflow DSL2, R/Bioconductor, HPC con Singularity + Slurm. Prioriza claridad, reproducibilidad y mantenibilidad por encima de completitud. Nada de sobre-ingeniería.

## Comunicación

- Chat en **español**, directo y técnico. Se asume background de bioinformática; no expliques lo básico.
- Comentarios de código en **inglés**, concisos, solo cuando la intención no es obvia. Tono de desarrollador experimentado, sin frases de IA.

  **Mal:** `# This function reads the input and returns the result`
  **Bien:** `# Spread multi-assigned reads 1/k across their miRNAs`

## Fuente de verdad

- `plan_idea.md` manda: alcance (§1), hechos (§2), decisiones (§3), parámetros (§5), estructura (§6), módulos (§8), salidas (§9), hitos (§10).
- **Antes de implementar un módulo, lee su sección §8** y respeta sus entradas/salidas exactas.
- Si el plan y el código no cuadran, pregunta en el chat en vez de improvisar.

## Alcance de cada orden

Una orden = **un módulo**. Toca solo los ficheros que se te indiquen. No adelantes fases, no refactorices módulos vecinos, no "aproveches" para arreglar cosas fuera de alcance. `main.nf` debe poder leerse de arriba abajo.

## Reglas duras

1. **No tests.** Sin `tests/`, `test/`, `bin/tests/` ni CI de tests. La validación es una **ejecución real** sobre el batch (plan §10).
2. **No crear carpetas ni ficheros fuera del árbol de abajo.** Si crees que hace falta uno, pregunta antes.
3. **No añadir** `nextflow_schema.json`, `subworkflows/`, `conf/` partido ni plantilla nf-core. Este pipeline es DSL2 plano por decisión (D2).
4. **Ningún umbral o ruta hardcodeado** en el código. Todo parámetro sale de `nextflow.config` y se vuelca a `run_info/pipeline_info/` (plan §5, D9).
5. Los datos de nf-core son **de solo lectura**. Todo se escribe en `outdir`.
6. **Nada específico del batch** en el código. Las diferencias entre muestras entran solo por `metadata.csv`, `contrasts.csv` y los flags de QC.
7. No crear ficheros `.md` nuevos salvo petición expresa. Las decisiones largas van a `docs/internal/`.
8. Actualiza `task.todo` (una fila) al cerrar cada orden.
9. Comentarios en inglés. Sin emojis.

## Árbol del repositorio

Este es **todo** el árbol permitido (plan §6). Marca lo que ya existe:

```
smrna-nextflow/
├── README.md                       # existe
├── AGENTS.md                       # este fichero
├── plan_idea.md                    # spec, existe
├── main.nf                         # por crear — flujo completo, legible de arriba abajo
├── nextflow.config                 # por crear — params, imagen, perfiles (singularity, docker, slurm)
├── modules/
│   ├── inputs.nf                   # M0, M1
│   ├── analysis.nf                 # M2, M3, M4
│   ├── reads.nf                    # M6 (fase 2)
│   └── extras.nf                   # M5, M7, M8 (fase 3)
├── bin/                            # un script de R por etapa, con --help (optparse)
│   ├── validate_inputs.R           # M0
│   ├── build_counts.R              # M1
│   ├── qc_composition.R            # M2
│   ├── deseq2_group.R              # M3
│   ├── summarize_de.R              # M4
│   ├── biotype_full.R              # M6
│   ├── isomir_profile.R            # M7
│   ├── functional.R                # M8
│   └── report.Rmd                  # M5
├── containers/smrna/Dockerfile     # existe — R/Bioconductor + bowtie + bedtools
├── .github/workflows/build-containers.yml   # existe — construye y publica en GHCR
├── assets/                         # plantillas genéricas: metadata.csv, contrasts.csv, params.yaml
├── demo/batch_2026_04/             # ejemplo real: metadata, contrasts, params.yaml, make_metadata.R, run.sh
└── docs/
    └── internal/                   # existe — papers, state_art.md, plan_idea_v1_original.md
```

Salidas de `nextflow run` (`results/`, `work/`, `.nextflow/`) van fuera del control de versiones. El `outdir` (`results/`) se organiza en `deliverables/` (tablas y figuras para el wet lab), `raw/` (volcados nativos de auditoría) y `run_info/` (parámetros, umbrales y `pipeline_info/`); el detalle está en `plan_idea.md` §9.

## Convenciones

**Nextflow**
- DSL2. Procesos en `UPPER_SNAKE_CASE`, canales en `ch_camelCase`.
- Todo proceso lleva: `tag`, `label` de recursos, `container params.container`, `publishDir` definido junto al proceso, `when`, `input`, `output`, `script`.
- Labels de recursos: `process_single`, `process_medium`, `process_high`.
- **Nunca `file()` en un `.config`**: no es DSL de Nextflow ahí. Usa `new File(...)`.
- Sin rutas absolutas en el código. Las referencias grandes se pasan por `params` (D8).
- Mensajes de error con qué falló, dónde, qué se esperaba y cómo arreglarlo.

**R (`bin/`)**
- `optparse` con `--help`; valida argumentos antes de ejecutar. Usa `message()`, `warning()`, `stop()`.
- `suppressPackageStartupMessages()` al cargar librerías.
- `set.seed(params.seed)` (D9).
- Figuras siempre en PDF y PNG.
- Tablas TSV: cabecera, tabuladores, sin comillas, punto decimal. Primera columna `sample_id` o `contrast_id` cuando aplique.

**Invocación de `bin/`**
- Llama los scripts como `Rscript ${projectDir}/bin/<script>.R ...` (patrón probado en `bulkAnnex`; requiere `singularity.autoMounts = true`). No los copies a `work/`.

## Validación (Definition of Done)

No hay tests. Cada orden termina con:

1. Verificación barata: `nextflow lint` (si está disponible) y/o `nextflow run . -stub-run` para comprobar el cableado.
2. Ejecución real:
   ```bash
   nextflow run . -params-file demo/batch_2026_04/params.yaml -profile singularity,slurm -resume
   ```
3. Pega en el chat las cabeceras de las tablas generadas, el recuento de filas y las rutas publicadas.

No se pasa a la fase siguiente sin haber ejecutado la anterior sobre el batch real.

## Contenedor

- `containers/smrna/Dockerfile` es la fuente de verdad. `nextflow.config` asigna una sola imagen por defecto a todos los procesos, sustituible con `params.container` (D3).
- Para actualizarla: editar el Dockerfile y hacer push a `main` (o `workflow_dispatch`); la Action construye y publica `ghcr.io/<owner>/smrna-nextflow`. El paquete debe ser público para que Singularity lo descargue.
- No metas dependencias nuevas sin avisar en el chat (qué hace falta y por qué).

## Estado del proyecto

- `task.todo` es la checklist viva: una fila por orden, con estado. Sirve para saber qué toca ahora.
- `plan_idea.md` es el spec congelado. No lo dupliques en `task.todo`.
- Commits pequeños: uno por módulo aceptado.
