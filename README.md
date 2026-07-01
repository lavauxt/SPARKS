<p align="center">
  <img src="assets/logo.png" alt="SPARKS logo" width="350"/>
</p>

<h1 align="center">SPARKS</h1>

<p align="center">
  Single-cell RNA-seq analysis toolkit
</p>

<p align="center">
  <a href="https://www.gnu.org/licenses/gpl-3.0">
    <img src="https://img.shields.io/badge/License-GPLv3-blue.svg" alt="GPL-3 License"/>
  </a>
</p>

**SPARKS** (Single-cell Pipeline for Analysis of RNA-seq Systems) is a modular R package designed to streamline scRNA-seq analysis. It automates the workflow from raw count matrices (Alevin or Cell Ranger) to cluster identification, annotation, sex scoring, and differential expression analysis using **Seurat**, **SingleR**, and **scDblFinder**.

It provides an end-to-end framework for preprocessing, quality control, dimensionality reduction, clustering, visualization, and biological interpretation of single-cell transcriptomics datasets.

Species currently supported: **Human** and **Mouse**.

## Key Features

- **Automated Workflow**: Processes multiple samples and comparison groups in a single run.
- **Flexible Input**: Supports both 10X Genomics (Cell Ranger) and Alevin output formats.
- **Comprehensive QC**: Automated filtering based on mitochondrial content, feature counts, and doublet detection.
- **Cell Type Annotation**: Integrated SingleR support with customizable reference datasets.
- **Pathway Enrichment**: Built-in support for GSEA/ssGSEA/UCell analysis via MSigDB.
- **Advanced Scoring**: Built-in modules for sex scoring and cell cycle regression.
- **Subset Analysis**: Targeted re-clustering and analysis for specific cell populations defined in the config.

## Installation

You can install the development version of SPARKS from GitHub:

```r
# Install devtools if not already installed
if (!require("devtools", quietly = TRUE)) install.packages("devtools")

# Install SPARKS (dependencies = TRUE also pulls in the optional report
# packages: rmarkdown, knitr, DT, plotly, htmltools, RColorBrewer)
devtools::install_github("lavauxt/SPARKS", dependencies = TRUE)
```

> **HTML reports (`QC_report_*.html` / `Results_report_*.html`)** are
> generated automatically at the end of each comparison group's run.
> They require the `rmarkdown` package (plus `DT`/`plotly`/`htmltools`/
> `RColorBrewer` for the interactive Results report). If you installed
> SPARKS without `dependencies = TRUE`, install them once with:
> ```r
> install.packages(c("rmarkdown", "knitr", "DT", "plotly", "htmltools", "RColorBrewer"))
> ```
> The `.Rmd` templates ship inside the installed package (`inst/rmd/`) and
> are located automatically — you no longer need to copy `qc_report.Rmd` /
> `results_report.Rmd` next to your config file for reports to generate.

## Quick Start

### 1. Prepare your Sample Table
Create a `sample_table.tsv` defining your data folders and experimental conditions:

| folder_id | protocol | comparison_group |
| :-------- | :------- | :--------------- |
| Sample_A  | WT       | Group_1          |
| Sample_B  | KO       | Group_1          |

### 2. Configure the Pipeline
Copy the provided mouse template and modify it for your study. You can override specific parameters (e.g., QC thresholds, resolution) in a separate YAML file.

### 3. Run the Pipeline

```{r run-sparks, eval=FALSE}
library(SPARKS)

# Run the complete analysis
# For mouse
sparks(
  base_config_path     = "./inst/config_template_mouse.yaml",
  override_config_path = "./inst/examples/config.yaml"
)

# For human
sparks(
  base_config_path     = "./inst/config_template_human.yaml",
  override_config_path = "./inst/examples/config.yaml"
)
```

## Configuration

The entire analytical pipeline is governed via a structured YAML configuration 
matrix. Below is a detailed mapping of the block structure parameters.

* Project Ingestion Landscape (pipeline)
Maps workspace directory contexts, target organisms, and metadata registry links.
```
  pipeline:
    data_dir: "./datas"             # Source directory tracking raw count roots
    results_dir: "./results"        # Target directory for all output files
    species_target: "Mouse"         # Organism directive ("Mouse" or "Human")
    sample_table: "./sample_table.tsv"  # File path to sample tracking matrix
```
*  Organism Profile Criteria (species)
Tracks sequence element clearance filters and establishes reference parameters 
for cell type taxonomy annotations.
```
  species:
    gene_removal_pattern: "^(Rps|Rpl|Rrn|Rn|Hb|Gm).*|.*Rik$" # Clears noise transcripts
    mt_pattern: "^mt-"                                      # Identifies mitochondrial genes
    genes_to_remove: ["Xist", "Tsix", "Kdm6a"]              # Explicit blacklist array
    ref_primary: "ImmGenData"                               # Primary SingleR reference
    ref_secondary: "MouseRNAseqData"                        # Secondary SingleR reference
```
*  Covariate Evaluation Systems (sex_scoring & cell_cycle)
Enables biometric state evaluations and directs downstream regression scaling.
```
  sex_scoring:
    run: true               # Calculates individual sex representation scores
    regress: true           # Regresses out sex-linked expression variation during scaling
    markers:
      female: ["Xist", "Tsix", "Kdm6a"]
      male: ["Eif2s3y", "Ddx3y", "Uty"]

  cell_cycle:
    run: true               # Quantifies cell cycle metrics
    regress: true           # Regresses out cell cycle signatures to avoid cell cycle confounding
```
*  Technical Quality Control Boundaries (qc)
Establishes rigid thresholds for filtering out low-quality cell events and 
transcript artifacts.
```
  qc:
    min_features: 200       # Minimum unique gene count required per cell
    max_counts: 30000       # Maximum acceptable depth profile tracking UMIs
    max_mt_percent: 10      # Percentage ceiling for mitochondrial transcript ratios
    min_cells: 3            # Baseline gene feature filtering requirement
```
*  Seurat Dimensional Scaling Matrix (processing)
Directs down-stream cluster processing, community partitioning, and dimensionality reductions.
```
  processing:
    vars_to_regress: ["percent.mt"]  # Non-genomic metadata covariates to regress out
    pca_dims_from: 1                 # Evaluation start index for PCA space
    pca_dims_to: 20                  # SNN graph build dimension ceiling
    npcs: 50                         # Total principal components computed initially
    cluster_resolution: 0.5          # Louvain modularity clustering resolution scale
    reduction: "umap"                # Low-dimensional embedding selection
```
*  Differential Gene Discovery (deg)
Configures strict statistical criteria used for identifying differentially expressed markers.
```
  deg:
    logfc_threshold: 0.25            # Minimum natural log fold-change cut-off
    min_pct: 0.1                     # Minimum feature detection frequency in target groups
    min_cells_per_group: 10          # Processing representation limit to qualify for testing
    min_p_val_adj: 0.01              # Adjusted significance ceiling limit
    avg_expression_layers: ["data", "counts"] # Layer extractions checked during analysis
```
*  Reference Classification Horizons (singler)
Maps classification targets, dictionary scopes, and matching conditions.
```
  singler:
    unassigned_label: "Unassigned"   # Label applied when correlation metrics match poorly
    min_cells_per_group: 10          # Population limits guiding assignment testing
    labels:
      - name: "singleR_labels_main"
        ref_field: "label.main"      # Broad lineage designation category mapping
        is_fine: false
      - name: "singleR_labels_fine"
        ref_field: "label.fine"      # High-resolution lineage designation mapping
        is_fine: true
```
*  Functional Enrichment Architecture (escape)
Orchestrates functional enrichment scoring routines via internal MSigDB parsing.
```
  escape:
    run: true
    method: "AUCell"                # Analytical engine selection ("AUCell", "UCell", "ssGSEA", "GSVA")
    library: ["H", "C5"]            # Focus databases targeted (e.g., Hallmark, Gene Ontology)
    min_size: 5                     # Minimum dimensional gene pathway footprint filter
```
*  Feature Target Registries (genes)
Directs automated expression tracing visualization outputs and coordinate matrix calculations.
```
  genes:
    genes_to_plot:                  # Global registry driving FeaturePlot/Violin generation
      - "Sphk1"
      - "Cd19"
      - "Cd3e"
    corr_genes_x: ["S1pr1", "S1pr2"] # Coordinate factors driving expression correlation matrices
    corr_genes_y: ["Ccl19", "Cxcl13"]
```
*  Directs automated expression tracing visualization outputs and coordinate matrix calculations.
```
  genes:
    genes_to_plot:                  # Global registry driving FeaturePlot/Violin generation
      - "Sphk1"
      - "Cd19"
      - "Cd3e"
    corr_genes_x: ["S1pr1", "S1pr2"] # Coordinate factors driving expression correlation matrices
    corr_genes_y: ["Ccl19", "Cxcl13"]
```
*  Gene Signature Panels (gene_signatures)
Named gene lists that get a dedicated per-cell z-score heatmap, an averaged
z-score heatmap, and a DotPlot (average expression + percent expressed)
automatically generated for **every grouping column** — `seurat_clusters`,
every SingleR label column, and every subset's own cluster/type column —
of **both** the Main analysis unit and every subset/subcluster defined below.
No extra wiring per subset is needed: add a signature once and it is applied
everywhere.
```
  gene_signatures:
    - name: "Lipid_Scavenger"
      genes: ["Cd36", "Ackr3", "Ackr4", "Plpp3", "Stab1", "Stab2", "Cav1", "Cav2", "Dab2"]
    - name: "Vascular_Tightness"
      genes: ["Cdh5", "Cldn5", "Pdgfa", "Pdgfb", "Tek", "Tjp1", "Wnt2", "Emcn"]
```
Output (per grouping column `<group_col>`, per analysis unit `<Main|SubsetName>`):
```
<Main|SubsetName>/<group_col>/Heatmap/
  ├── <Unit>_<group_col>_<Signature>_Heatmap_PerCell.png    # z-score per cell
  ├── <Unit>_<group_col>_<Signature>_Heatmap_Aggregated.png # z-score averaged per group
  └── <Unit>_<group_col>_<Signature>_DotPlot.png            # avg expression + % expressed
```
Genes not present in the dataset are skipped automatically (with a message);
a signature is skipped entirely if fewer than 2 of its genes are found, or if
the grouping column has fewer than 2 groups.

*  Declarative Sub-Clustering Directives (subsets)
Defines cell population routing rules for automated secondary extraction, 
re-clustering, and sub-population profiling.
```
  subsets:
    - display_name: "Fibroblasts"
      match_col: "singleR_labels_main"   # Targeting metadata column
      pattern: ["Fibroblasts"]           # Target pattern used for routing extraction
      exact_match: true                  # Strict array text validation matching
      type_col: "fibroblast_type"        # Assigned sub-clustering metadata tracking slot
      pca_dims_from: 1                   # Dedicated sub-space calculation dimensions start
      pca_dims_to: 10                    # Dedicated sub-space calculation dimensions stop
      color: "darkgreen"                 # Color mapping tracking parameter
      genes: ["Ccl19", "Cxcl13"]         # Subpopulation targeted gene expression tracking set
      label_rules: []                    # Optional rule modifiers
```
## Output Structure

The pipeline generates an organized results directory:
```
results_dir/
  └── Group_1/
      ├── QC/
      │   ├── VlnPlot_Group_1.png            # Ingestion profile distributions
      │   ├── ScatterPlot_Group_1.png        # Depth vs feature count tracking
      │   └── QC_report_Group_1.html         # Auto-generated QC HTML report
      ├── Main/                              # (repeated per subset, e.g. Fibroblasts/, Endothelial/)
      │   └── seurat_clusters/               # (repeated per grouping: singleR labels, type_col, ...)
      │       ├── UMAP/
      │       │   ├── UMAP_Main_seurat_clusters.png
      │       │   └── UMAP_DEG_Main_seurat_clusters.png
      │       ├── DEG/
      │       │   └── AllMarkers_Main_seurat_clusters.txt
      │       ├── VlnPlot/
      │       └── Heatmap/
      │           ├── ZScore_TopGenes_Main_seurat_clusters.png
      │           ├── ExprHeatmap_Main_seurat_clusters.png
      │           ├── Main_seurat_clusters_Lipid_Scavenger_Heatmap_PerCell.png
      │           ├── Main_seurat_clusters_Lipid_Scavenger_Heatmap_Aggregated.png
      │           ├── Main_seurat_clusters_Lipid_Scavenger_DotPlot.png
      │           ├── Main_seurat_clusters_Vascular_Tightness_Heatmap_PerCell.png
      │           ├── Main_seurat_clusters_Vascular_Tightness_Heatmap_Aggregated.png
      │           └── Main_seurat_clusters_Vascular_Tightness_DotPlot.png
      ├── Escape/
      │   ├── Enrichment_Scores.tsv          # Pathway scores per single cell
      │   └── Pathway_Heatmap.png            # Cluster-averaged pathway enrichment profiles
      ├── Results_report_Group_1.html        # Auto-generated interactive Results HTML report
      └── RData/
          └── Group_1_Seurat_Processed.rds   # Production-ready Seurat binary object file
```
## Changelog

- **Fixed**: HTML reports (`QC_report_*.html`, `Results_report_*.html`) were
  silently never generated once SPARKS was installed as a package, because
  the `.Rmd` templates were not shipped inside the package and the lookup
  only checked the config's directory and the current working directory.
  Templates now live in `inst/rmd/` and are located via `system.file()` as a
  final fallback. `rmarkdown`, `knitr`, `DT`, `plotly`, `htmltools`, and
  `RColorBrewer` are now declared in `Suggests` so `dependencies = TRUE`
  installs them.
- **Fixed**: `DESCRIPTION`'s `Collate` field referenced `pipeline_core.R` and
  `escape_analysis.R`, which don't exist in this package — the real files are
  `main.R` and `enrichment_analysis.R`. Corrected to avoid collation issues.
- **Added**: `gene_signatures` config block — named gene panels (e.g. a lipid
  scavenger receptor set, a vascular tightness/junction set) that
  automatically get a per-cell z-score heatmap, an averaged z-score heatmap,
  and a DotPlot for every grouping column of the Main unit and every subset.

## License

This project is licensed under the [GNU General Public License v3.0](LICENSE).