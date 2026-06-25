#!/usr/bin/env python3
"""09_scrna.py — localize the 11 shared hub genes to cell types in heart &
skeletal muscle using CZ CELLxGENE Census. Streams X with axis_query
(out-of-core), accumulating per-cell-type sum/nonzero counts for the 11 genes —
avoids materializing a huge AnnData and the slow obs_coords path.
Outputs per tissue:
  results/09_scrna/09_hub_celltype_<tissue>.csv
  results/09_scrna/dotplot_<tissue>.png
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "09_scrna"); os.makedirs(OUT, exist_ok=True)
LOG = os.path.join(ROOT, "logs", "run.log")
def log(m):
    import datetime
    s = f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] 09 | {m}"
    print(s, flush=True); open(LOG, "a").write(s + "\n")

import cellxgene_census, tiledbsoma as soma

HUB = ["EIF4EBP1","PITPNM1","MAGED2","NUDT4","MTFP1","PPM1K","NAMPT",
       "CCND1","FBLN1","ERAP2","CENPV"]
CENSUS_VERSION = "2023-12-15"
TISSUES = {
    "heart":  "tissue_general == 'heart'",
    "muscle": "tissue_general == 'musculature' or tissue == 'skeletal muscle tissue'",
}

def run_tissue(census, name, tfilter):
    hs = census["census_data"]["homo_sapiens"]
    base = f"({tfilter}) and is_primary_data == True"
    var = (census["census_data"]["homo_sapiens"]["ms"]["RNA"]["var"]
           .read(value_filter="feature_name in [%s]" % ",".join(f"'{g}'" for g in HUB))
           .concat().to_pandas())
    vmap = dict(zip(var["soma_joinid"], var["feature_name"]))  # var joinid -> gene
    q = hs.axis_query("RNA",
            obs_query=soma.AxisQuery(value_filter=base),
            var_query=soma.AxisQuery(value_filter="feature_name in [%s]"
                                     % ",".join(f"'{g}'" for g in HUB)))
    obs = q.obs(column_names=["soma_joinid", "cell_type"]).concat().to_pandas()
    n = len(obs)
    if n == 0:
        log(f"{name}: 0 cells — skipping"); return
    vc = obs["cell_type"].value_counts()
    keep = list(vc[vc >= 200].index)
    n_by_ct = vc.to_dict()
    # integer codes for vectorized accumulation
    ct_code = {c: i for i, c in enumerate(keep)}            # cell types we keep
    obs["code"] = obs["cell_type"].map(ct_code)             # NaN for dropped types
    code_ser = pd.Series(obs["code"].values, index=obs["soma_joinid"].values)
    gene_list = list(HUB)
    g_code = {gj: gene_list.index(vmap[gj]) for gj in vmap}  # var joinid -> gene index
    g_ser = pd.Series(g_code, dtype="float64")              # for vectorized reindex
    nC, nG = len(keep), len(gene_list)
    smat = np.zeros((nC, nG)); zmat = np.zeros((nC, nG))
    log(f"{name}: {n:,} cells, {vc.size} cell types ({nC} with >=200); streaming X (vectorized) for {nG} hubs")

    for batch in q.X("raw").tables():
        oj = batch["soma_dim_0"].to_numpy(); vj = batch["soma_dim_1"].to_numpy()
        dv = np.log1p(batch["soma_data"].to_numpy())
        cc = code_ser.reindex(oj).to_numpy()               # cell-type code per entry (NaN if dropped)
        gc = g_ser.reindex(vj).to_numpy()                  # gene index per entry (NaN if not a hub)
        m = ~np.isnan(cc) & ~np.isnan(gc)
        cc = cc[m].astype(int); gc = gc[m].astype(int); dv = dv[m]
        np.add.at(smat, (cc, gc), dv)
        np.add.at(zmat, (cc, gc), 1.0)
    q.close()

    rows = []
    for c in keep:
        i = ct_code[c]; ncell = n_by_ct[c]
        for j, g in enumerate(gene_list):
            rows.append(dict(tissue=name, cell_type=c, gene=g,
                             mean_log1p=smat[i, j] / ncell,
                             pct_expr=100.0 * zmat[i, j] / ncell,
                             n_cells=int(ncell)))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, f"09_hub_celltype_{name}.csv"), index=False)
    log(f"{name}: wrote table ({res.cell_type.nunique()} cell types x {len(HUB)} hubs)")

    # dotplot (matplotlib): x=gene, y=cell type, size=pct, color=mean_log1p
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        # order cell types by total hub expression; cap to top 25 for readability
        order = (res.groupby("cell_type")["mean_log1p"].sum().sort_values(ascending=False))
        cts_keep = list(order.head(25).index)
        d = res[res.cell_type.isin(cts_keep)].copy()
        d["cell_type"] = pd.Categorical(d.cell_type, categories=cts_keep[::-1], ordered=True)
        d["gene"] = pd.Categorical(d.gene, categories=HUB, ordered=True)
        fig, ax = plt.subplots(figsize=(1.1 + 0.55 * len(HUB), 1 + 0.34 * len(cts_keep)))
        sca = ax.scatter(d.gene.cat.codes, d.cell_type.cat.codes,
                         s=4 + d.pct_expr * 2.2, c=d.mean_log1p, cmap="Reds",
                         edgecolors="grey", linewidths=.3)
        ax.set_xticks(range(len(HUB))); ax.set_xticklabels(HUB, rotation=55, ha="right", fontsize=8)
        ax.set_yticks(range(len(cts_keep))); ax.set_yticklabels(cts_keep[::-1], fontsize=7)
        ax.set_title(f"{name}: shared-hub expression by cell type (CELLxGENE Census)", fontsize=10)
        cb = fig.colorbar(sca, ax=ax, shrink=.5); cb.set_label("mean log1p", fontsize=8)
        ax.margins(.04); fig.tight_layout()
        fig.savefig(os.path.join(OUT, f"dotplot_{name}.png"), dpi=130, bbox_inches="tight")
        log(f"{name}: dotplot saved ({len(cts_keep)} cell types shown)")
    except Exception as e:
        log(f"{name}: dotplot failed: {e}")

def main():
    log(f"scRNA hub localization via Census {CENSUS_VERSION} (streaming); {len(HUB)} hubs")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        for name, tf in TISSUES.items():
            try:
                run_tissue(census, name, tf)
            except Exception as e:
                import traceback; log(f"{name}: FAILED {type(e).__name__}: {e}")
                traceback.print_exc()
    log("09 done.")

if __name__ == "__main__":
    main()
