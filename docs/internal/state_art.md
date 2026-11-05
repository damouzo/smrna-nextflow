¡Hola! Es un gran salto pasar del *bulk RNA-seq* clásico a la bioinformática de *small RNAs*. Mientras que en bulk RNA te mueves en el terreno de los transcriptos completos (mRNAs y lncRNAs), aquí entramos en un universo microscópico (habitualmente de 15 a 50 nucleótidos) donde las reglas cambian por completo: la estructura tridimensional, las modificaciones químicas y la fragmentación de RNAs estructurales (como tRNAs y rRNAs) mandan tanto o más que la simple complementariedad de secuencias.

A continuación, te presento una guía detallada del estado del arte en este campo, conectando los conceptos con el output de tu pipeline **`nf-core/smrnaseq`** (con datos de Novogene) y las últimas revoluciones científicas reflejadas en la literatura actual.

---

## 1. El ABC de un análisis estándar de Small RNA (vs. Bulk RNA)

Si vienes de bulk RNA, muchos conceptos te sonarán familiares, pero la implementación técnica en *small RNA* tiene particularidades críticas:

* **QC (Control de Calidad) y Adapter Trimming:** En *small RNA*, los reads de secuenciación son a menudo más cortos que el propio adaptador de Illumina. Por tanto, el **adapter trimming** no es opcional, es el paso más crítico: si el adaptador no se elimina perfectamente, las lecturas no alinearán contra nada. Herramientas como *Trim Galore* (integradas en `nf-core/smrnaseq`) se encargan de esto.
* **Alineamiento y Anotación (El gran cambio frente a HISAT2 + featureCounts):**
* *En bulk RNA:* Aligns al genoma $\rightarrow$ cuantificación por genes (ej. `featureCounts`).
* *En *small RNA* (Small RNA-seq):* Un enfoque puramente genómico falla porque hay cientos de copias idénticas de tRNAs o rRNAs dispersas por el genoma. Los pipelines especializados usan **estrategias de alineamiento jerárquico** (usualmente con mapas tipo Bowtie o STAR). Primero se busca emparejar contra bases de datos de referencia específicas (miRBase para miRNAs, GtRNAdb para tRNAs, Rfam, etc.) antes de mapear al genoma de referencia.




* **Normalización:** Olvídate del TPM o FPKM basado en la longitud del gen (aquí casi todos los small RNAs tienen longitudes muy similares). Se suele utilizar **RPM (Reads Per Million)** normalizados al total de reads mapeados (o al total de reads de *non-coding RNA*), o bien factores de normalización robustos (como TMM de *edgeR*) si cuentas con suficientes réplicas biológicas.
* **DGE (Differential Gene Expression / Abundance):** A nivel estadístico, herramientas como *DESeq2* o *edgeR* se siguen utilizando para encontrar especies diferencialmente abundantes. Sin embargo, debes tener cuidado con el **sesgo de composición de la librería**: si una clase masiva de small RNAs (como los *tsRNAs*) cambia drásticamente entre condiciones, puede sesgar la normalización global.
* **Análisis Funcional / GSEA:** En bulk RNA buscas enriquecimiento en rutas KEGG o Gene Ontology (GO). En *small RNA*, el análisis funcional es muy distinto:
* Para *miRNAs*, realizas predicción de dianas (*target prediction* con herramientas como TargetScan o miRDB) y luego enriquecimiento de sus mRNAs blancos.
* Para *non-canonical small RNAs* (como tsRNAs o snoRNAs), el análisis funcional a menudo implica buscar **interacciones con proteínas (RBPs)**, cambios en modificaciones post-transcripcionales o análisis de estructuras secundarias.





---

## 2. El Universo de las Especies de Small RNA

Históricamente, el análisis de *small RNA* se cegaba al 100% en los **miRNAs** (microRNAs, ~21-23 nt, encargados del silenciamiento post-transcripcional mediante el complejo RISC/AGO). Hoy sabemos que los miRNAs representan a menudo **menos del 5%** del total de small RNAs celulares. El resto pertenece a un ecosistema fascinante:

### A. tsRNAs / tDRs (tRNA-derived small RNAs)

Son fragmentos generados a partir de los tRNAs (ya sea de los extremos 5' o 3', o mitades completas de tRNA).

* *Función:* Lejos de ser "basura" o degradación celular, actúan bajo estrés celular (hipoxia, daño, etc.). Por ejemplo, en el reciente trabajo de **Li et al. (2025)**, se descubre que un tDR específico del extremo 3' (tRNA-Asp-GTC-3'tDR) se dispara durante la hipoxia renal, forma estructuras de tipo **G-quadruplex (G4)**, secuestra a la enzima modificadora de RNA **PUS7** y desencadena un proceso de degradación selectiva llamado *RNA autophagy* para proteger al riñón.



### B. snoRNAs y sdRNAs (Small Nucleolar RNAs & snoRNA-derived small RNAs)

* *Función canónica:* Tradicionalmente confinados al nucléolo, los snoRNAs (clases C/D y H/ACA) guían modificaciones químicas esenciales en los rRNAs (2'-O-metilación y pseudouridilación) y participan en la biogenesis de ribosomas.


* *Funciones no canónicas:* Como detalla la exhaustiva revisión de **Wang et al. (2025)**, los snoRNAs también regulan el procesamiento del extremo 3' de los mRNAs, modulan el *splicing* alternativo e interaccionan con lncRNAs. Además, se fragmentan en piezas más pequeñas llamadas **sdRNAs**, los cuales actúan de manera muy similar a los miRNAs regulando la traducción o la estabilidad de mRNAs.



### C. rsRNAs (rRNA-derived small RNAs)

Fragmentos derivados de los rRNAs ribosomales. Durante años los bioinformáticos los filtraban automáticamente tachándolos de "contaminación de laboratorio" o ruido ribosomal. Hoy se sabe que muchos son estables, altamente regulados y cumplen funciones de señalización.

### D. piRNAs (PIWI-interacting RNAs)

De mayor tamaño (21-35 nt), altamente enriquecidos en células germinales y gónadas, esenciales para reprimir la actividad de los retrotransposons.

### E. Otros non-canonical sncRNAs

Fragmentos procedentes de *Y RNAs* (ysRNAs), *vault RNAs* (vtsRNAs), *snRNAs* (snsRNAs) e incluso fragmentos de ARMs mensajeros (*msRNAs*).

---

## 3. Las Nuevas Fronteras: ¿Qué está apareciendo en la Bioinformática de Small RNA?

Si vas a procesar datos de Novogene con pipelines estándar (`smrnaseq`), debes conocer los dos grandes cuellos de botella tecnológicos que la bioinformática actual está resolviendo:

1. **El problema de las modificaciones químicas y los extremos bloqueados:**
Los *small RNAs* no clásicos (tsRNAs, rsRNAs) están plagados de modificaciones post-transcripcionales (como metilaciones $m^1A$, $m^3C$, $m^1G$) y extremos no convencionales (como fosfatos 2',3'-cíclicos o extremos 5'-OH). Las librerías de *small RNA* tradicionales de Novogene **son ciegas ante ellos**, ya que la transcriptasa inversa se bloquea al tropezar con una metilación y los adaptores no ligan bien si los extremos no son 5'-P / 3'-OH.


* *La solución actual:* Protocolos avanzados como **PANDORA-seq** (que incorporan tratamientos enzimáticos previos con T4 PNK y la demetilasa AlkB) combinados con pipelines bioinformáticos específicos como **SPORTS1.1** permiten destapar este "litoral oscuro" del transcriptoma.




2. **Resolución a nivel de secuencia individual (Sequence-level resolution):**
Herramientas modernas como *SPORTS1.1* o *seqpac* no se limitan a darte un conteo por "gen", sino que analizan la **coordenada exacta de inicio y fin** del fragmento sobre el tRNA o rRNA parental. Esto es vital, porque un mismo tRNA puede generar un fragmento protector (en el extremo 3') o un fragmento tóxico/pro-apoptótico (en el extremo 5') dependiendo de la patología.


3. **Biología Estructural computacional aplicada:**
Como se vio en los descubrimientos sobre tsRNAs, la función de muchos small RNAs no depende de la complementariedad de bases lineal, sino de su capacidad de plegamiento tridimensional (ej. la formación de **G-quadruplexes / G4**). Integrar predicciones estructurales (con paquetes como *NUPACK*) en el flujo bioinformático es una tendencia al alza.



---

## 4. ¿Cómo encaja todo esto con tu output de `smrnaseq` y Novogene?

Tu ejecución de `nf-core/smrnaseq` (v2.2.4) con datos de Novogene te generará una estructura típica de Nextflow:

* Informes de **MultiQC** con métricas de calidad y tasas de asignación.
* Tablas de conteos cuantificadas principalmente para **miRNAs** (consultando bases de datos como miRBase).
* Alineamientos mapeados contra el genoma de referencia.

**Una advertencia importante para tu transición desde Bulk RNA:**
Si tus datos de Novogene provienen de una librería comercial estándar de *small RNA* sin protocolos de retirada de modificaciones (como PANDORA-seq), la inmensa mayoría de tus lecturas que mapeen a tRNAs o rRNAs pueden ser artefactos de truncamiento de la transcriptasa inversa o estar subrepresentadas. Si tu objetivo biológico principal son los miRNAs, el pipeline estándar te servirá de maravilla. Pero si investigas estrés celular, patologías complejas o respuestas tisulares (como el eje de tsRNAs y autofagia descrito en el riñón por *Li et al.*), ten en cuenta que los métodos estándar solo rascan la superficie de lo que tus datos de Novogene podrían albergar.

¿Te gustaría profundizar en alguno de estos puntos, por ejemplo, en cómo interpretar los archivos de salida de anotación o cómo configurar un pipeline específico para la búsqueda de tsRNAs en tus datos?