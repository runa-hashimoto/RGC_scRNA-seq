from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse



# 1. Settings and paths

H5AD_PATH = Path(
    "./1_scRNAseq_preprocessing_and_RGC_extraction_results/"
    "subset_cluster3_for_tricycle.h5ad"
)

GENESET_PATH = Path(
    "./2_QuickGO_gene_set_generation_results/"
    "merged_union_gene_list_mouse.tsv"
)

OUTDIR = Path("./3_filter_GO_genes_by_RGC_expression_results")
OUTDIR.mkdir(parents=True, exist_ok=True)

MIN_EXPR_PCT = 0.05
EXPR_THRESHOLD = 0



# 2. Read input data

if not H5AD_PATH.exists():
    raise FileNotFoundError(f"H5AD file not found: {H5AD_PATH}")

if not GENESET_PATH.exists():
    raise FileNotFoundError(f"Gene set file not found: {GENESET_PATH}")

adata_rgc = sc.read_h5ad(H5AD_PATH)

if adata_rgc.n_obs == 0:
    raise ValueError("The RGC dataset contains no cells.")

gene_set_df = pd.read_csv(GENESET_PATH, sep="\t")

if "gene_symbol" not in gene_set_df.columns:
    raise ValueError(
        f"'gene_symbol' column not found in gene set file. "
        f"Columns: {list(gene_set_df.columns)}"
    )

go_genes = (
    gene_set_df["gene_symbol"]
    .dropna()
    .astype(str)
    .unique()
    .tolist()
)



# 3. Match GO genes to the RGC dataset

genes_in_adata = [
    gene for gene in go_genes
    if gene in adata_rgc.var_names
]

genes_not_in_adata = [
    gene for gene in go_genes
    if gene not in adata_rgc.var_names
]



# 4. Calculate expression frequency in RGCs

X = adata_rgc[:, genes_in_adata].X

if sparse.issparse(X):
    X = X.toarray()
else:
    X = np.asarray(X)

expr_bool = X > EXPR_THRESHOLD

n_cells_expressed = expr_bool.sum(axis=0)
pct_cells_expressed = n_cells_expressed / adata_rgc.n_obs
mean_expr_all_cells = X.mean(axis=0)

mean_expr_positive_cells = []

for i in range(X.shape[1]):
    vals = X[:, i]
    vals_pos = vals[vals > EXPR_THRESHOLD]

    mean_expr_positive_cells.append(
        vals_pos.mean() if len(vals_pos) > 0 else 0
    )



# 5. Select genes expressed in at least 5% of RGCs

result_df = pd.DataFrame(
    {
        "gene_symbol": genes_in_adata,
        "n_cells_expressed": n_cells_expressed,
        "pct_cells_expressed": pct_cells_expressed,
        "mean_expr_all_RGC_cells": mean_expr_all_cells,
        "mean_expr_positive_RGC_cells": mean_expr_positive_cells,
    }
)

result_df["expressed_in_RGC"] = (
    result_df["pct_cells_expressed"] >= MIN_EXPR_PCT
)

result_df = result_df.sort_values(
    by=[
        "expressed_in_RGC",
        "pct_cells_expressed",
        "n_cells_expressed",
        "mean_expr_all_RGC_cells",
    ],
    ascending=[False, False, False, False],
).reset_index(drop=True)

expressed_df = result_df[
    result_df["expressed_in_RGC"]
].copy()

not_expressed_df = result_df[
    ~result_df["expressed_in_RGC"]
].copy()



# 6. Save results

result_df.to_csv(
    OUTDIR / "GO_union_genes_expression_summary_cluster3_RGC_all.csv",
    index=False,
)

expressed_df.to_csv(
    OUTDIR / "GO_union_genes_expressed_in_cluster3_RGC_pct05.csv",
    index=False,
)

not_expressed_df.to_csv(
    OUTDIR / "GO_union_genes_not_expressed_in_cluster3_RGC_pct05.csv",
    index=False,
)

expressed_df[["gene_symbol"]].to_csv(
    OUTDIR / "GO_union_genes_expressed_in_cluster3_RGC_pct05_gene_symbols_only.csv",
    index=False,
)

pd.DataFrame(
    {"gene_symbol": genes_not_in_adata}
).to_csv(
    OUTDIR / "GO_union_genes_not_found_in_adata.csv",
    index=False,
)



# 7. Summary

print("===== Summary =====")
print("GO union genes:", len(go_genes))
print("Genes found in RGC dataset:", len(genes_in_adata))
print("Genes not found in RGC dataset:", len(genes_not_in_adata))
print("RGC cells:", adata_rgc.n_obs)
print("Expression threshold:", EXPR_THRESHOLD)
print("Minimum RGC expression fraction:", MIN_EXPR_PCT)
print("Genes expressed in >=5% of RGCs:", expressed_df.shape[0])
print("Genes expressed in <5% of RGCs:", not_expressed_df.shape[0])
print("Output directory:", OUTDIR.resolve())