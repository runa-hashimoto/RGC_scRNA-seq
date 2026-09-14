### E14.5 mouse cortical RGC scRNA-seq analysis
This repository contains the code used to reanalyze publicly available scRNA-seq data from E14.5 mouse cortex, with a focus on radial glial cells (RGCs).
#
### Workflow
Run the scripts in numerical order:
1. `1_scRNAseq_preprocessing_and_RGC_extraction.py`
2. `2_QuickGO_gene_set_generation.R`
3. `3_filter_GO_genes_by_RGC_expression.py`
4. `4_cell_cycle_trajectory_analysis.R`


### Input data
The scRNA-seq dataset is available from GEO:`GSM4635075`

Set the input file path in `1_scRNAseq_preprocessing_and_RGC_extraction.py`.


### Gene set used for the manuscript analysis
The fixed 883-gene set used for the manuscript analysis is provided in:
`GO_union_genes_expressed_in_cluster3_RGC_pct05_gene_symbols_only.csv`

Because QuickGO is continuously updated, rerunning the GO retrieval may produce a slightly different gene set.


### Environment
Python environment information is provided in:
`environment_scenv.yml`
