from pathlib import Path
import os
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
from matplotlib.colors import LinearSegmentedColormap



# 0. Settings and paths

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor="white")
scvi.settings.seed = 0


H5_PATH = Path(
    "/path/to/GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5"
)

OUTDIR = Path("./1_scRNAseq_preprocessing_and_RGC_extraction_results")
FIGDIR = OUTDIR / "figures_pdf"


OUTDIR.mkdir(parents=True, exist_ok=True)
FIGDIR.mkdir(parents=True, exist_ok=True)

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["savefig.bbox"] = "tight"

BATCH_LABEL = "E14.5"
CLUSTER_KEY = "leiden_scVI_0.4"
SOLO_CUTOFF = 0.80

if not H5_PATH.exists():
    raise FileNotFoundError(
        f"Input file not found: {H5_PATH.resolve()}"
    )

print("Output directory:", OUTDIR.resolve())
print("Figure directory:", FIGDIR.resolve())



# 1. Read 10x HDF5 and preserve raw counts

adata = sc.read_10x_h5(H5_PATH, gex_only=True)
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()

adata.obs["batch"] = BATCH_LABEL
adata.obs["batch"] = adata.obs["batch"].astype("category")

print(f"Raw: cells={adata.n_obs:,}, genes={adata.n_vars:,}")



# 2. Cell and gene filtering

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=20)

adata.var["mt"] = adata.var_names.str.match(r"(?i)^mt-")
sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, max_counts=22000)
sc.pp.filter_cells(adata, min_genes=500)
adata = adata[adata.obs["pct_counts_mt"] < 5].copy()

print(f"After QC: cells={adata.n_obs:,}, genes={adata.n_vars:,}")
if adata.n_obs != 3267 or adata.n_vars != 12367:
    warnings.warn(
        "QC output differs from the original run "
        f"(expected 3267 cells x 12367 genes; observed {adata.n_obs} x {adata.n_vars})."
    )



# 3. Normalize and log-transform for expression visualization

sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
sc.pp.log1p(adata)
adata.raw = adata



# 4. scVI latent representation from raw count data

scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch",
)
model = scvi.model.SCVI(adata)
model.train()
adata.obsm["X_scVI"] = model.get_latent_representation()



# 5. SOLO doublet detection and cutoff 0.80 filtering

solo = scvi.external.SOLO.from_scvi_model(model)
solo.train()
pred_soft = solo.predict(soft=True)

pred_soft.index = pred_soft.index.str.replace(r"-0$", "", regex=True)

if "doublet" not in pred_soft.columns:
    raise ValueError(
        "SOLO soft prediction does not contain a 'doublet' column. "
        f"Columns: {list(pred_soft.columns)}"
    )

common_cells = adata.obs_names.intersection(pred_soft.index)
adata = adata[common_cells].copy()
adata.obs["SOLO_doublet_score"] = (
    pred_soft.reindex(adata.obs_names)["doublet"].astype(float)
)
adata.obs["SOLO_prediction_0.80"] = np.where(
    adata.obs["SOLO_doublet_score"] >= SOLO_CUTOFF,
    "doublet",
    "singlet",
).astype(str)

adata = adata[adata.obs["SOLO_prediction_0.80"] == "singlet"].copy()

print(f"Cells after SOLO cutoff 0.80: {adata.n_obs:,}")
if adata.n_obs != 3005:
    warnings.warn(
        "The number of cells after SOLO cutoff 0.80 differs from the original run "
        f"(expected 3005; observed {adata.n_obs})."
    )



# 6. Neighborhood graph, UMAP, and Leiden clustering

sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_scVI")
sc.tl.umap(adata, min_dist=0.5)
sc.tl.leiden(
    adata,
    key_added=CLUSTER_KEY,
    resolution=0.4,
)

adata.write(OUTDIR / "E14_5_after_SOLO_cutoff0.80_Leiden0.4.h5ad")


cluster_key = "leiden_scVI_0.4"
n_top_genes = 1000

marker_outdir = OUTDIR / "leiden_scVI_0.4_top1000_marker_genes"
marker_outdir.mkdir(exist_ok=True)

adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

sc.tl.rank_genes_groups(
    adata,
    groupby=cluster_key,
    method="wilcoxon",
    use_raw=False,
    key_added="rank_genes_" + cluster_key
)

all_marker_list = []

for cluster in adata.obs[cluster_key].cat.categories:

    df = sc.get.rank_genes_groups_df(
        adata,
        group=cluster,
        key="rank_genes_" + cluster_key
    )

    df_top = df.head(n_top_genes).copy()
    df_top.insert(0, "cluster", cluster)
    df_top.insert(
        1,
        "rank",
        range(1, df_top.shape[0] + 1)
    )

    all_marker_list.append(df_top)

    df_top.to_csv(
        marker_outdir /
        f"cluster_{cluster}_top1000_marker_genes.csv",
        index=False
    )

all_markers_df = pd.concat(
    all_marker_list,
    axis=0,
    ignore_index=True
)

all_markers_df.to_csv(
    marker_outdir /
    "leiden_scVI_0.4_all_clusters_top1000_marker_genes.csv",
    index=False
)

print(all_markers_df.head())



# 7. Figure: Leiden-cluster UMAP

sc.pl.umap(
    adata,
    color=CLUSTER_KEY,
    frameon=False,
    show=False,
)
plt.savefig(FIGDIR / "Fig_scRNA_Leiden0.4_UMAP.pdf", format="pdf", bbox_inches="tight")
plt.close()



# 8. Figure: marker-gene dot plot

marker_dict = {
    "RGC": ["Pax6", "Sox2", "Vim", "Hes1", "Hes5", "Fabp7", "Slc1a3", "Nes"],
    "IP": ["Eomes", "Neurod1", "Neurod4", "Insm1", "Btg2"],
    "Proliferating": ["Mki67", "Top2a", "Cenpf", "Ube2c"],
    "EN": ["Neurod2", "Tbr1", "Satb2", "Dcx", "Tubb3", "Rbfox3",
        "Bcl11b", "Cux1", "Cux2", "Rnd2", "Neurod6"],
    "IN": ["Gad1", "Gad2", "Dlx1", "Dlx2", "Lhx6", "Arx", "Maf", "Mafb"],
    "CR": ["Reln", "Calb2", "Cxcl12"],
    "E/V": ["Pecam1", "Rgs5", "Pdgfrb"],
}

marker_dict_present = {
    group: [gene for gene in genes if gene in adata.var_names]
    for group, genes in marker_dict.items()
}
marker_dict_present = {
    group: genes for group, genes in marker_dict_present.items() if genes
}

print("Marker genes used in the dot plot:")
for group, genes in marker_dict_present.items():
    print(group, genes)

var_names_flat = []
var_group_positions = []
var_group_labels = []
start = 0
for group_name, genes in marker_dict_present.items():
    end = start + len(genes) - 1
    var_names_flat.extend(genes)
    var_group_positions.append((start, end))
    var_group_labels.append(group_name)
    start = end + 1

dp = sc.pl.DotPlot(
    adata,
    var_names=var_names_flat,
    groupby=CLUSTER_KEY,
    standard_scale="var",
    var_group_positions=var_group_positions,
    var_group_labels=var_group_labels,
    figsize=(max(14, len(var_names_flat) * 0.32), 6),
)
dp.make_figure()
ax_dict = dp.get_axes()
main_ax = ax_dict["mainplot_ax"]
main_ax.set_xticklabels(
    main_ax.get_xticklabels(),
    rotation=90,
    ha="center",
    va="top",
)
dp.fig.subplots_adjust(top=0.82)
dp.fig.savefig(
    FIGDIR / "Fig_scRNA_marker_dotplot.pdf",
    format="pdf",
    bbox_inches="tight",
)
plt.close(dp.fig)



# 9. Figure: marker-gene UMAPs

gray_to_black = LinearSegmentedColormap.from_list(
    "gray_to_black",
    ["#d9d9d9", "#000000"],
)

marker_dict = {
    "RGC": [
        "Pax6", "Sox2", "Vim", "Hes1", "Hes5",
        "Fabp7", "Slc1a3", "Nes"
    ],
    "IPC": [
        "Eomes", "Neurod1", "Neurod4", "Insm1", "Btg2"
    ],
    "CellCycle": [
        "Mki67", "Top2a", "Cenpf", "Ube2c"
    ],
    "Excitatory neuron": [
        "Neurod2", "Tbr1", "Satb2", "Dcx", "Tubb3",
        "Rbfox3", "Bcl11b", "Cux1", "Cux2",
        "Neurod1", "Rnd2", "Neurod6"
    ],
    "Inhibitory neuron": [
        "Gad1", "Gad2", "Dlx1", "Dlx2",
        "Lhx6", "Arx", "Maf", "Mafb"
    ],
    "CR": [
        "Reln", "Calb2", "Cxcl12"
    ],
    "Endothelial / vascular": [
        "Pecam1", "Rgs5", "Pdgfrb"
    ]
}

marker_all = []
for genes in marker_dict.values():
    marker_all.extend(genes)

marker_all = list(dict.fromkeys(marker_all))

marker_all_present = [
    gene for gene in marker_all
    if gene in adata.var_names
]

print(
    f"Marker genes detected: "
    f"{len(marker_all_present)} / {len(marker_all)}"
)
print(marker_all_present)

sc.pl.umap(
    adata,
    color=marker_all_present,
    frameon=False,
    cmap=gray_to_black,
    ncols=4,
    show=False
)

plt.savefig(
    FIGDIR / "Fig_scRNA_all_marker_genes_UMAP.pdf",
    format="pdf",
    bbox_inches="tight"
)

plt.close()



# 10. Cell-type annotation used in the paper

cluster_to_celltype = {
    "0": "EN",
    "1": "EN",
    "4": "EN",
    "2": "IP",
    "6": "IP",
    "3": "RGC",
    "5": "IN",
    "8": "CR",
    "7": "E/V",
}

celltype_order = ["RGC", "IP", "EN", "IN", "CR", "E/V"]
adata.obs["celltype_annot"] = adata.obs[CLUSTER_KEY].astype(str).map(cluster_to_celltype)
adata.obs["celltype_annot"] = pd.Categorical(
    adata.obs["celltype_annot"],
    categories=celltype_order,
    ordered=True,
)

unassigned = sorted(
    adata.obs.loc[adata.obs["celltype_annot"].isna(), CLUSTER_KEY]
    .astype(str)
    .unique()
)
if unassigned:
    raise ValueError(f"Unassigned Leiden clusters: {unassigned}")

print("Cell numbers by annotation:")
print(adata.obs["celltype_annot"].value_counts(dropna=False))
print("Cluster x annotation:")
print(pd.crosstab(adata.obs[CLUSTER_KEY], adata.obs["celltype_annot"], dropna=False))



# 11. Figure: annotated cell-type UMAP

sc.pl.umap(
    adata,
    color="celltype_annot",
    legend_loc="on data",
    frameon=False,
    show=False,
)
plt.savefig(FIGDIR / "Fig_scRNA_celltype_UMAP.pdf", format="pdf", bbox_inches="tight")
plt.close()



# 12. Extract RGC dataset for downstream tricycle analysis

cluster_key = "leiden_scVI_0.4"

if cluster_key not in adata.obs.columns:
    raise ValueError(f"{cluster_key} not found in adata.obs")


adata_sub = adata[
    adata.obs[cluster_key].astype(str).isin(["3"])
].copy()

print("RGC cells:", adata_sub.n_obs)
print(adata_sub.obs[cluster_key].value_counts())


if "X_scVI" not in adata_sub.obsm:
    raise ValueError("X_scVI not found in the RGC dataset.")

sc.pp.neighbors(
    adata_sub,
    n_neighbors=10,
    use_rep="X_scVI",
)

sc.tl.umap(
    adata_sub,
    min_dist=0.5,
)


sc.pl.umap(
    adata_sub,
    color=cluster_key,
    frameon=False,
    show=False,
)

plt.savefig(
    FIGDIR / "Fig_scRNA_RGC_recomputed_UMAP.pdf",
    format="pdf",
    bbox_inches="tight",
)

plt.close()


adata_sub.write(
    OUTDIR / "subset_cluster3_for_tricycle.h5ad"
)

print(
    "RGC dataset for tricycle saved:",
    OUTDIR / "subset_cluster3_for_tricycle.h5ad",
)



# 13. Serotonin-receptor subtype screening and violin plots

serotonin_receptors = [
    "Htr1a", "Htr1b", "Htr1d", "Htr1f",
    "Htr2a", "Htr2b", "Htr2c",
    "Htr3a", "Htr3b", "Htr3c", "Htr3d", "Htr3e",
    "Htr4", "Htr5a", "Htr5b", "Htr6", "Htr7",
]

receptor_screen = pd.DataFrame(
    {
        "gene": serotonin_receptors,
        "retained_after_filtering": [gene in adata.var_names for gene in serotonin_receptors],
    }
)
receptor_screen.to_csv(OUTDIR / "serotonin_receptor_screening.csv", index=False)
print("Serotonin receptor screening:")
print(receptor_screen)

receptors_present = [gene for gene in serotonin_receptors if gene in adata.var_names]
print("Serotonin receptor genes available for downstream visualization:")
print(receptors_present)

clusters = adata.obs[CLUSTER_KEY].astype(str)
plot_group_order = [
    "Radial glial cell",
    "Intermediate progenitor",
    "Cajal-Retzius cell",
    "Excitatory neuron",
    "Inhibitory neuron",
]
plot_group_to_cells = {
    "Radial glial cell": adata.obs_names[clusters.isin(["3"])],
    "Intermediate progenitor": adata.obs_names[clusters.isin(["2", "6"])],
    "Cajal-Retzius cell": adata.obs_names[clusters.isin(["8"])],
    "Excitatory neuron": adata.obs_names[clusters.isin(["0", "1", "4"])],
    "Inhibitory neuron": adata.obs_names[clusters.isin(["5"])],
}


def get_gene_expression_vector(adata_obj, gene_name):
    x = adata_obj[:, gene_name].X
    if sparse.issparse(x):
        return np.asarray(x.toarray()).ravel()
    return np.asarray(x).ravel()


sns.set(style="whitegrid", context="notebook")

for gene in receptors_present:
    expr = get_gene_expression_vector(adata, gene)
    gene_df = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "expr": expr,
        }
    )

    rows = []
    for group_name in plot_group_order:
        cell_ids = plot_group_to_cells[group_name]
        tmp = gene_df.loc[gene_df["cell_id"].isin(cell_ids)].copy()
        tmp["group"] = group_name
        rows.append(tmp)

    plot_df = pd.concat(rows, axis=0, ignore_index=True)
    plot_df["group"] = pd.Categorical(
        plot_df["group"],
        categories=plot_group_order,
        ordered=True,
    )

    fig, ax = plt.subplots(figsize=(7.5, 5))

    sns.violinplot(
        data=plot_df,
        x="group",
        y="expr",
        order=plot_group_order,
        inner=None,
        cut=0,
        ax=ax,
    )

    sns.stripplot(
        data=plot_df,
        x="group",
        y="expr",
        order=plot_group_order,
        jitter=0.20,
        size=3,
        alpha=0.45,
        ax=ax,
    )

    ax.set_xticklabels(plot_group_order, rotation=55, ha="right")
    ax.set_title(gene, fontstyle="italic")
    ax.set_xlabel("")
    ax.set_ylabel("Expression (a.u.)")
    plt.tight_layout()

    fig.savefig(
        FIGDIR / f"Fig_serotonin_receptor_{gene}_violin.pdf",
        format="pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

print(f"Analysis completed. Results: {OUTDIR.resolve()}")

