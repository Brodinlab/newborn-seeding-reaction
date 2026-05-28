# newborn-seeding-reaction
Code to reproduce figures from 'A seeding reaction to colonizing microbes in human newborns'

### File structure

Scripts for generating figures and tables are kept in `scripts/`.<br>
Input and output files sorted by respective datatype

```
├── scripts
│   ├── newborn-seeding-reaction.ipynb
│   └── ..R
├── input
│   ├── cytof
│   ├── IgGseq
│   ├── luminex
│   ├── metadata
│   ├── mass-spec
│   ├── metagenomics
│   ├── nCounter
│   ├── olink
│   ├── nulisa
│   ├── olink_functional
│   ├── olink_nulisa
│   └── scatac
│   └── rnaseq
└── output
    ├── cytof
    ├── IgGseq
    ├── luminex
    ├── metadata
    ├── mass-spec
    ├── metagenomics
    ├── nCounter
    ├── olink
    ├── nulisa
    ├── olink_functional
    ├── olink_nulisa
    └── scatac
    └── rnaseq
```    


All analyses were performed using R 4.4.2


#### `newborn-seeding-reaction.ipynb`
Notebook with the complete code for reproducing figures with metagenomic, olink, nulisa, olink_functional, mass-spec, luminex, IgGeq, cytof, nCounter and metadata.<br>


#### Bulk and scRNA-seq figure scripts

| Script | Input | Output |
|--------|-------|--------|
| `scripts/bulk_mrnaseq.r` | `input/bulkRNAseq/` | `output/bulkRNAseq/` (Fig4D, Fig5D, Fig5E) |
| `scripts/scRNAseq.r` | `input/scRNAseq/` + external UMI (see `DATA_SOURCES.md`) | `output/scRNAseq/` (Fig7C–7E) |

```bash
Rscript scripts/bulk_mrnaseq.r

# scRNA: re-plot if you have a local cache (gitignored); otherwise use output/scRNAseq/*.pdf
Rscript scripts/scRNAseq.r

# scRNA: build local cache from UMI (~30–60 min)
REBUILD_FROM_RAW=TRUE Rscript scripts/scRNAseq.r
```

`output/scRNAseq/cache/` is gitignored (author-local RDS after a full rebuild). The repo ships `output/scRNAseq/*.pdf`; cloning does not require cache.

### Relevant R packages
```
ggpubr 0.6.0
ArchR 1.0.3
DESeq2 1.48.1
clusterProfiler 4.16.0
limma 3.64.1
EnhancedVolcano 1.16.0
Tidyverse 2.0.0
readxl 1.4.5
kml 2.5.0
ggplot2 3.5.1
ggrepel 0.9.6
microshades 1.13
dplyr 1.1.4
tidyr 1.3.1
reshape2 1.4.4
stringr 1.5.1
phyloseq 1.50.0
mia 1.14.0
gee 4.13.29
arrow 18.1.0.1
factoextra 1.0.7
ALDEx2 1.38.0
ggtree 3.14.0
tidytree 0.4.6
ggtreeExtra 1.16.0
ggstar 1.0.4
treeio 1.30.0
vegan 2.7.0
ComplexHeatmap 2.22.0
circlize 0.4.16
OlinkAnalyze 4.2.0
ggpp 0.5.8.1
CellGrid 0.6.1 
FlowSOM 2.2.0
vite 0.4.10
ggraph 2.1.0
```
