# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a biohackathon project (KIDS23-Team7) for spatial transcriptomics analysis, focusing on integrating and analyzing Akoya PhenoCycler data using Seurat v5. The project includes generic analysis functions, integration pipelines, and a Shiny web interface for interactive data exploration.

## Technology Stack

- **Primary Language**: R
- **Key Packages**:
  - Seurat v5 (with v5 assay objects enabled via `options(Seurat.object.assay.version = "v5")`)
  - data.table for data I/O
  - harmony for batch integration
  - shiny/shinyWidgets for web interface
  - openxlsx for cell type annotation tables
- **Spatial Platforms Supported**: Akoya PhenoCycler, Visium, Vizgen MERSCOPE, Nanostring CosMx

## Architecture

### Core Function Library: GenericFunctions.R

This file contains the main analysis pipeline with reusable functions:

1. **Data Loading**: `csv_to_seurat()` - Converts expression matrices and metadata CSVs to Seurat objects with spatial FOV (field of view) coordinates
2. **Normalization & Sketching**: `norm_and_sketch()` - Handles platform-specific normalization (CLR for Akoya, SCTransform for others) and downsampling via SketchData
3. **Analysis Pipeline**: `seurat_workflow()` - Complete workflow including PCA/UMAP, neighbor finding, and clustering with multiple resolutions
4. **Cell Type Annotation**: `annotate_cell_types()` - Uses sc-type algorithm to annotate clusters based on marker gene sets from Excel files

### Integration Scripts

- **Group1_integration_cellAnno.R** & **Group2_integration_cellAnno.R**: Sample integration workflows using Harmony for batch correction across multiple samples
- **TestSample_Akoya_seurat5.R**: Single-sample processing template for Akoya data with hexagonal coordinate transformation

### Shiny Application: test_ybae/app.R

Interactive web interface for uploading data, running QC, and performing analysis. Sources GenericFunctions.R.

## Common Development Tasks

### Running R Scripts

This is an R project with RStudio configuration. Open in RStudio or run R scripts directly:

```r
# Source the generic functions first
source("GenericFunctions.R")

# Then use in your scripts or interactively
obj <- csv_to_seurat("path/to/expression.csv", "path/to/metadata.csv")
```

### Running the Shiny App

```r
# From the project root
shiny::runApp("test_ybae/app.R")
```

## Key Implementation Patterns

### Seurat v5 Assay Objects

Always set this option at the start of scripts:
```r
options(Seurat.object.assay.version = "v5")
```

### Spatial Coordinate Handling

Akoya data requires hexagonal grid coordinate transformation:
```r
r <- 1/2
R <- (2 / sqrt(3)) * r
ak_meta$y <- ak_meta$Y - min(ak_meta$Y) + 1
ak_meta$x <- ak_meta$X - min(ak_meta$X) + 1
ak_meta$y <- ak_meta$y * R * (3/2)
ak_meta$x <- (ak_meta$x + 1) / 2
ak_meta$y <- -ak_meta$y
```

### Multi-Sample Integration Workflow

Standard pattern used in integration scripts:
1. Load individual Seurat objects
2. Set DefaultAssay to 'Akoya'
3. Add Sample metadata
4. Merge objects
5. Re-sketch the merged data
6. Run Harmony integration on 'Sample' variable
7. Use 'harmony' reduction for downstream analysis (UMAP, clustering)

### Cell Type Annotation with sc-type

External functions are sourced from GitHub:
```r
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
```

Cell type tables are Excel files with columns: `tissueType`, `cellName`, `geneSymbolmore1` (comma-separated gene lists)

## Data Assumptions

- Expression matrices: Rows = cells, Columns = markers/genes
- Metadata: Must include Object_ID (cell IDs), X, Y coordinates for spatial data
- Cell type reference tables: Excel format with gene symbol lists
- Default assay names: 'Akoya', 'Spatial', 'Vizgen', 'Nanostring', 'Xenium'
- Sketched assay is always named 'sketch'

## Analysis Parameters

- Default PCA dimensions: 30
- Default clustering resolutions: 0.3, 0.5, 1, 1.5, 2
- Default sketch size: 50,000 cells (via LeverageScore method)
- Features treated as all variable when panel size < 1000
- sc-type confidence threshold: scores < ncells/5 marked as "Unknown"
