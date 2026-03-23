# GeNetME

**Gene-network meta-analysis of prioritised genes from ME/CFS studies**

GeNetME is an R workflow that merges independently published ME/CFS gene modules into a unified protein–protein interaction (PPI) network, computes network metrics, and performs over-representation analysis (ORA) across four pathway and ontology databases. It is part of the [MetaME](https://github.com/paolomaccallini-hub/MetaME) project.

---

## Requirements

| Software | Version tested |
|----------|----------------|
| R | 4.4.1 (2024-06-14 ucrt) |
| RStudio | 2026.1.1.403 |
| OS | Windows 11 24H2 (build 10.0.26200.8037) |

### R packages

**CRAN**

```r
install.packages(c(
  "data.table",
  "ggplot2",
  "GOplot",
  "httr",
  "igraph",
  "jsonlite",
  "magick",
  "Matrix",
  "readxl",
  "stringr"
))
```

**Bioconductor**

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "clusterProfiler",
  "DOSE",
  "enrichplot",
  "org.Hs.eg.db",
  "pathview",
  "ReactomePA",
  "TissueEnrich"
))
```

---

## Installation

Clone the repository and set the working directory to the repository root.

```bash
git clone https://github.com/paolomaccallini-hub/GeNetME.git
```

In RStudio: Session → Set Working Directory → To Source File Location.

No compilation or additional setup is required.

---

## External data

The following files are downloaded and cached automatically on first run. An internet connection is required for the first execution only.

| File | Source | Size (approx.) | Local path |
|------|--------|----------------|------------|
| `9606.protein.links.v12.0.txt.gz` | STRING v12.0 | 400 MB | `Data/STRING/` |
| `9606.protein.info.v12.0.txt.gz` | STRING v12.0 | 5 MB | `Data/STRING/` |
| `gene_info.gz` | NCBI Gene | 1 GB (filtered to human and saved as `human_gene_info.gz`) | `Data/NCBI/` |
| `media-9.xlsx` | Zhang S. et al. 2025 (medrxiv) | — | `Data/` |
| `media-2.xlsx` | Sardell JM. et al. 2025 (medrxiv) | — | `Data/` |

A background universe file (`Data/My_Universe.tsv`) mapping all STRING genes to NCBI Entrez IDs is also built on first run. This step is slow (≈ 20,000 lookups against the local NCBI database) but runs only once.

---

## Usage

The workflow is split into two files.

| File | Role |
|------|------|
| `GeNetME_Func.R` | Library loading, global parameters, data download and caching, all function definitions |
| `GeNetME_Main.R` | Analysis script — sources `GeNetME_Func.R` and runs the full pipeline |

Run the full pipeline by sourcing the main script:

```r
source("GeNetME_Main.R")
```

The pipeline writes a checkpoint to `PF_output/Disease_results_list.rds` (disease module + adjacency matrix). The `saveRDS` / `readRDS` block is left intentionally explicit so the analysis can be resumed from that point without repeating the network-construction and STRING API steps.

### Key parameter

`STRING.co` in `GeNetME_Func.R` sets the STRING combined-score threshold used for both network construction and the ORA background universe:

```r
STRING.co <- 0.4  # medium confidence (default)
```

---

## Input gene modules

| Module | Cases | Sequencing | Gene-mapping method | Criteria | Genes | Reference |
|--------|-------|------------|---------------------|----------|-------|-----------|
| Zhang S. 2025 | 464 | WGS | Deep learning on rare variants | ICC-IOM | 115 | [PMC12047926](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047926/) |
| Sardell JM. 2025 | 14,767 | Axiom UKB array | Combinatorial analysis of DecodeME GWAS | CCC-IOM | 259 | [medrxiv 2025.12.01.25341362](https://www.medrxiv.org/content/10.64898/2025.12.01.25341362v2) |

Zhang module: genes filtered at `q_value < 0.02`. Sardell module: all genes from Supplementary Table 3, sheet 3. Both modules are further filtered to genes present in STRING v12.0 before network construction.

---

## Methods

1. Gene names are mapped to STRING preferred names via the STRING `/get_string_ids` API.
2. NCBI Entrez Gene IDs are assigned from a local copy of NCBI `gene_info`, querying canonical symbols first, then synonyms.
3. Genes appearing in both modules are collapsed to a single row; `list.name` records both source labels and `list.count` records the number of contributing modules.
4. A symmetric weighted adjacency matrix is built from STRING v12.0 PPI data filtered at `STRING.co`. Symmetry is verified explicitly before returning.
5. Graph metrics (connected components, degree, strength, normalised score) are computed with igraph.
6. ORA is performed against four databases. Background universe: all STRING genes with a valid Entrez ID (~19,338). Multiple testing correction: Benjamini–Hochberg; significance threshold: adjusted p ≤ 0.05.

| R package | Function | Database |
|-----------|----------|----------|
| clusterProfiler | `enrichKEGG()` | KEGG |
| clusterProfiler | `enrichGO()` | Gene Ontology (Cellular Component) |
| ReactomePA | `enrichPathway()` | Reactome |
| DOSE | `enrichDO()` | Disease Ontology |
| TissueEnrich | `teEnrichment()` | Human Protein Atlas |

7. KEGG pathway maps are rendered with pathview and moved to `PF_output/ORA/KEGG/`.
8. Reactome pathway gene lists are exported as plain text files for submission to the [Reactome PathwayBrowser](https://reactome.org/PathwayBrowser/#TOOL=AT).

---

## Output

```
Output/
├── Disease_results_list.rds          # Checkpoint: disease module data frame + adjacency matrix
├── All_genes_graph.tiff / .jpeg      # PPI network figure
├── All_genes_cytoscape.tsv           # Edge list for Cytoscape import
└── ORA/
    ├── ORA_Disease_Module.tsv        # Full ORA results table (all databases)
    ├── TORA_Disease_Module.tsv       # Tissue enrichment results table
    ├── Tissue_ORA.tiff / .jpeg       # Tissue enrichment bar chart
    ├── GO/
    │   ├── GO_ORA_goplot.tiff / .jpeg
    │   └── GO_ORA_cnetplot.tiff / .jpeg
    ├── KEGG/
    │   └── hsa*.Disease_ORA.png      # KEGG pathway maps (pathview, one per term)
    └── Reactome/
        └── R-HSA-*.txt               # Entrez ID lists for Reactome PathwayBrowser
```

`all_genes.tsv` (repository root) lists all 369 genes with their NCBI Entrez ID, STRING preferred name, source module(s), connected component, degree, strength, and normalised score.

---

## Repository files

| File | Description |
|------|-------------|
| `GeNetME_Func.R` | Functions and setup |
| `GeNetME_Main.R` | Main analysis script |
| `all_genes.tsv` | Disease module: 369 genes with network metrics |
| `All_genes_cytoscape.tsv` | Edge list for Cytoscape |
| `ORA_Disease_Module.tsv` | Full ORA results |
| `LICENSE` | GPL-3.0 |

---

## License

GPL-3.0. See [LICENSE](LICENSE).

---
