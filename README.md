### E14.5 mouse cortical RGC scRNA-seq analysis
This repository contains code for the reanalysis of publicly available scRNA-seq data from E14.5 mouse cortex, with a focus on radial glial cells (RGCs).

#

### Input data
The scRNA-seq data analyzed here were originally generated in a previous study and deposited in the Gene Expression Omnibus (GEO) under accession `GSM4635075`.


### Gene set used for the cell-cycle trajectory analysis
The 883-gene set used in the analysis reported in the paper is provided in this repository as `GO_union_genes_expressed_in_cluster3_RGC_pct05_gene_symbols_only_2026-04-24.csv`.

Because Gene Ontology annotations are updated over time, rerunning the QuickGO retrieval may produce a slightly different gene set.


### Environment
Python environment information is provided in `environment_scenv.yml`.

R environment information for the cell-cycle trajectory analysis is provided in `RGC_GO_pct05_sessionInfo.txt`.
