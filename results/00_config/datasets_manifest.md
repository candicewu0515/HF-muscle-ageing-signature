# Dataset manifest — HF × sarcopenia shared biomarker

Verified via NCBI E-utilities (db=gds) on 2026-06-25. Raw JSON:
`dataset_verification.json`, `dataset_verification_sar.json`.

## Decision: MAIN LINE holds (no plan B needed)

Both disease sides have a discovery cohort large enough for WGCNA
(rule of thumb: per-group N > 15–20) plus independent external-validation
cohorts. Sarcopenia is operationalized as **aging skeletal muscle (old vs
young)** — the field-standard proxy given bulk-data limits; flag explicitly
in Discussion as a limitation.

## HF side (heart, left ventricle)

| Role | Accession | N | Platform | Notes |
|------|-----------|---|----------|-------|
| Discovery (DEG + WGCNA) | **GSE57338** | 313 | GPL11532 (Affy Human Gene ST, array) | 177 HF LV + 136 non-failing LV. Large N → stable WGCNA. NOTE: brief said "RNA-Seq"; it is actually an **array** dataset. |
| External validation (ROC) | GSE5406 | 210 | GPL96 (Affy U133A) | ischemic + idiopathic CM vs non-failing |
| External validation (ROC) | GSE116250 | 64 | GPL16791 (RNA-seq) | HF human LV, RNA-seq → cross-platform check |
| Backup validation | GSE1145 / GSE16499 / GSE76701 | 107 / 30 / 8 | mixed | only if needed |
| scRNA (step 09) | GSE183852 | TBD | snRNA | human DCM/HF LV (Koenig 2022). **Re-verify before step 09** (NCBI DNS flaky now). |

## Sarcopenia side (skeletal muscle)

| Role | Accession | N | Platform | Notes |
|------|-----------|---|----------|-------|
| Discovery (DEG + WGCNA) | **GSE8479** | 65 | GPL2700 (Sentrix Human-6) | old vs young (+pre/post exercise). Baseline old≈young≈25 → per-group > 20. |
| External validation (ROC) | GSE1428 | 22 | GPL96 | titled "Skeletal muscle sarcopenia" (old sarcopenic vs young) |
| External validation (ROC) | GSE25941 | 36 | GPL570 | effect of age on muscle transcriptome |
| Backup / pooled WGCNA | GSE9103 / GSE38718 / GSE28392 / GSE47881 | 40/22/70/89 | all GPL570 | same platform → ComBat-poolable if more N needed |
| scRNA (step 09) | **GSE167186** | 93 | snRNA (GPL20301;24676) | aged muscle single-nuclei — resolves "muscle scRNA missing" risk |

## Red line — DO NOT USE

| Accession | Why |
|-----------|-----|
| GSE56815 | blood monocytes, pre/post-menopausal females (BMD/osteoporosis). Confirmed N=80, GPL96. Occupied by sarcopenia×osteoporosis niche. |

## Open items
- [ ] Re-verify GSE183852 (heart snRNA) before step 09 — NCBI DNS intermittently failing this session.
- [ ] Confirm GSE8479 exact group labels (young/old, pre/post) when parsing phenotype.
- [ ] Confirm GSE57338 etiology breakdown (ischemic vs DCM) in phenotype — pool as "HF" vs "control".
