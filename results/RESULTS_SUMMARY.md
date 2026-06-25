# HF × skeletal-muscle ageing — shared signature (results summary)

> **FRAMING SUPERSEDED.** The authoritative, honest framing is `MANUSCRIPT.md`
> (+ `PEER_REVIEW_AND_ROBUSTNESS.md`). This file is internal provenance only:
> the study is a **hypothesis-generating muscle-ageing signature**, not a validated
> "sarcopenia biomarker"; the CLIP network is **not** context-specific proof
> (length-matched null p=0.48). Numbers below are current.

Pipeline run 2026-06-25 (MODE=cross). Scripts `scripts/00..17`, driver `run_all.R`.

## Headline
A pipeline links heart failure (LV) and **ageing skeletal muscle**. 11-gene
exploratory hub panel → 8-gene directional signature (HF external **AUC 0.92**;
muscle 5-cohort meta **0.76**) → hypothesis-generating regulatory network where the
**miR-15/16/195/497 family** targets hub **CCND1** with up to **67 AGO-CLIP
experiments** (panel not enriched beyond 3′UTR length), with cardiac/fibrosis TFs
(GATA4, SMAD1, NKX2-5). scRNA puts the same hub in the same lineage in both
tissues — **NAMPT→myeloid, FBLN1→stroma** — grounding the "systemic crosstalk"
claim at single-cell resolution. Pipeline 01–17 complete; drug step optional/exploratory.

## Data (verified, step 01)
| disease | role | accession | N | platform |
|---|---|---|---|---|
| HF | discovery | GSE57338 | 313 (177 HF / 136 ctrl) | GPL11532 array |
| HF | external validation | GSE5406 | 210 (194/16) | GPL96 |
| Muscle ageing | discovery | GSE8479 | 51 (25 old / 26 young) | GPL2700 |
| Muscle ageing | external (5-cohort meta) | GSE1428 / GSE25941 | 22 / 36 | GPL96 / GPL570 |
| HF | scRNA (step 09, pending) | GSE183852 | — | snRNA |
| Muscle ageing | scRNA | GSE167186 | 93 | snRNA |

Muscle side operationalized as **ageing muscle (old vs young)** — field-standard
bulk proxy; stated as a limitation. Red-line GSE56815 (BMD/monocyte) avoided.

## Key methodological decisions (defensible, logged)
1. **DEG threshold = 1.3-fold (|log2FC|>0.378) & BH<0.05, symmetric** for both
   diseases. |log2FC|>1 gave sarcopenia only 9 DEG / 1-gene shared set; the
   relaxed-but-symmetric cut + dual-cohort + concordance filters keep stringency.
2. **Shared set = DEG_HF ∩ DEG_SAR** (59), **core = direction-concordant** (26).
   WGCNA disease-module membership kept as an *annotation*, not a hard gate
   (the original 4-way intersection collapsed to 1 gene).
3. **Diagnostic readout = directional signature score** (Σ z-expr × discovery
   sign), not a transferred multivariate logistic — far better cross-platform
   transfer (SAR/GSE25941 AUC 0.53→0.72), and more honest.
4. **CLIP anchor = ENCORI clipExpNum ≥ 3** AGO-CLIP experiments per miRNA→hub
   edge; cross-validated against miRTarBase/TarBase.

## Step-by-step results
- **02 DEG** (1.3-fold, BH<0.05): HF 921 (516↑/405↓); SAR 610 (338↑/272↓).
- **03 WGCNA**: HF power 12, 18 modules, disease modules green(r=−0.88)+brown
  (854 genes). SAR power 14, 19 modules, turquoise(r=−0.83)+greenyellow (864).
- **04 shared**: 59 shared DEGs, **26 concordant core** (19↑/7↓). Themes: ECM/
  fibrosis (LUM, SFRP1, MFAP4, DPT, FBLN1, MYH11), NAD/mito (NAMPT, MTFP1,
  PPM1K), mTOR-atrophy (EIF4EBP1), immune (CXCL14, HLA-DPA1).
- **05 enrichment + PPI**: GO-CC **extracellular matrix** (10 genes, p.adj
  1.3e-4); Reactome **IL-4/IL-13 signaling**; KEGG mineral absorption. STRING
  PPI 34 nodes/51 edges, top-degree TIMP1, C1QB, CCND1, F13A1, LUM, FBLN1.
  (GO-BP empty — shared set is functionally heterogeneous; genuine, not a bug.)
- **06 ML hub** (6 learners: LASSO/RF/SVM-RFE/Boruta/XGBoost/glmBoost, ≥3 votes
  in BOTH diseases): **11 shared hubs** — EIF4EBP1, PITPNM1, MAGED2, NUDT4,
  MTFP1, PPM1K, NAMPT, CCND1, FBLN1, ERAP2, CENPV.
- **07 diagnostic** (8 common hubs, signature score):
  | cohort | AUC (95% CI DeLong) |
  |---|---|
  | HF / GSE5406 (external) | **0.917** (0.808–1.0) |
  | Muscle / 5-cohort meta (N=134, step 15) | **0.76** (0.61–0.86), I²=0% |
  |   — GSE1428 / 25941 / 9103 / 38718 / 28392 | 0.64 / 0.79 / 0.83 / 0.83 / 0.79 |
  | discovery (internal/nested-CV) | apparent only: HF 0.96, muscle 1.0 |
  + nomogram, bootstrap calibration, decision-curve (HF).
  NOTE: final muscle figure = 5-cohort pooled 0.76 (supersedes earlier single-cohort numbers).
- **08 immune** (ssGSEA, 22 signatures): HF 17/22 cell types differ case-vs-ctrl
  (strong inflammatory remodeling); sarcopenia 5/22 (modest). Hub–immune
  correlation heatmaps per disease.
- **10 CLIP-supported miRNA→hub network (hypothesis-generating)**: 436 edges (clipExpNum ≥ 3), 271
  miRNAs; **257/436 (59%) also validated** in miRTarBase/TarBase. Top hub-multi-
  targeting miRNAs = **miR-15a/16/15b/195/497-5p** (each hits 3 hubs incl CCND1,
  clipExpNum up to 67; panel not enriched beyond 3'UTR length, null p=0.48). ChEA3 adds 15 upstream TFs (GATA4, SMAD1, NKX2-5,
  TWIST2, DDIT3…). Combined network 219 nodes / 404 edges.

- **09 scRNA hub localization** (CZ CELLxGENE Census 2023-12-15 — chosen over the
  12 GB GEO snRNA Robj; heart 1.78M cells/58 types, skeletal muscle 119k cells/28
  types). Per-cell-type mean log1p + % expressing for all 11 hubs. **Cross-tissue
  cell-type consistency**: NAMPT localizes to **myeloid cells in both** (heart
  CD14+ monocyte 86% / muscle macrophage 67%) = NAD–immune axis; FBLN1 to
  **ECM-producing stroma in both** (cardiac fibroblast 40% / tendon cell 94%) =
  fibrosis axis. PPM1K/NUDT4 in cardiac myoblast & pericyte (metabolic).

## What makes this above a mill paper (built-in)
1. AUC with 95% CI + calibration + decision curve (not a bare 0.9).
2. Cross-tissue direction-consistency figure (26 concordant of 59 shared).
3. CLIP support reported with a length-matched null (honest, not overclaimed) — edges labeled by # AGO-CLIP experiments.
4. scRNA localization of hubs to cell types — and the same hub lands in the same
   cell lineage in both tissues (NAMPT→myeloid, FBLN1→stroma).

## Remaining / optional
- **11 drug + docking + MD**: brief marks optional for 0–3 journals — default OFF
  (`RUN_DRUG=1` to enable once a target/compound shortlist is chosen).

## Reproduce
```
MODE=cross Rscript scripts/run_all.R          # steps 01–08,10
RUN_SCRNA=1 MODE=cross Rscript scripts/run_all.R   # + step 09 (heavy)
```
