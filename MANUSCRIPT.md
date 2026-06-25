# A shared cardiac–skeletal-muscle ageing transcriptional signature in heart failure implicates a miR-15/16–CCND1 regulatory axis


---

## Abstract

**Background.** Skeletal-muscle wasting frequently accompanies heart failure (HF)
and worsens prognosis, but whether the two tissues share a definable molecular
programme — and whether any such overlap is merely shared ageing — is unclear.

**Methods.** We integrated end-stage failing left-ventricular myocardium
(GSE57338, n=313) and ageing skeletal muscle (GSE8479, n=51) by independent
differential-expression and WGCNA analyses, intersected the disease signatures,
prioritised hub genes with six machine-learning selectors, and evaluated a
directional signature score in independent cohorts. Upstream regulation was mapped
with AGO-CLIP-supported miRNA interactions (ENCORI) and transcription factors
(ChEA3); hub genes were localised in 1.9 million single cells (CZ CELLxGENE
Census). Robustness was assessed by permutation, threshold-sensitivity,
age-adjustment and etiology-stratified analyses.

**Results.** Heart and muscle shared 59 dysregulated genes — significantly more
than chance (hypergeometric p=1.0×10⁻³) — enriched for extracellular-matrix,
NAD/immune-metabolic and mTOR programmes. An 8-gene signature discriminated HF in
an independent cohort (AUC 0.92, 95% CI 0.81–1.0) and, across five independent
ageing-muscle cohorts (N=134), gave a pooled AUC of 0.76 (95% CI 0.61–0.86,
I²=0%). The HF component was **robust to age** (AUC 0.95 after age adjustment;
signature p=1.8×10⁻¹⁸ vs age-alone AUC 0.64) and **etiology** (ischemic 0.96,
dilated 0.97; cross-etiology r=0.92), whereas the muscle signal weakened when
ageing genes were removed (non-ageing subset AUC 0.54–0.70), i.e. remained
ageing-linked. The predicted regulatory network was dominated by the canonical
**miR-15/16/195/497 → CCND1** cell-cycle axis (though not enriched beyond 3′UTR
length at the panel level); single-cell mapping placed NAMPT in myeloid cells and
FBLN1 in stromal cells in both tissues.

**Conclusions.** HF and skeletal-muscle ageing share a reproducible
tissue-remodelling and NAD/immune-metabolic transcriptional programme; the HF
component of this signature is robust to age and etiology, whereas the muscle
signal remains substantially ageing-linked. We provide this as a
hypothesis-generating map of a candidate cardio-muscle ageing axis, not a
validated diagnostic biomarker.

---

## Introduction

Heart failure (HF) is frequently complicated by loss of skeletal-muscle mass and
function — a continuum spanning sarcopenia and cardiac cachexia [3] — that
independently predicts mortality [1,2]. Mechanistic links proposed between failing
myocardium and wasting muscle include systemic inflammation (TNF-α, IL-6),
oxidative stress, mitochondrial dysfunction, activation of the
ubiquitin–proteasome and autophagy systems (FoxO–MuRF1/atrogin-1), insulin
resistance and RAAS activation. Whether these converging stresses produce a
*shared, definable transcriptional programme* across the two tissues — and whether
any such overlap reflects disease-specific cross-talk rather than the chronological
ageing common to both populations — has not been tested directly.

Two features motivate a cross-tissue design. First, the tissues are anatomically
distinct (myocardium vs skeletal muscle), so an overlap derived from independent
per-tissue analyses speaks to *systemic* rather than local biology. Second, the
miRNA layer (myomiRs miR-1/133/206/208/499; fibrosis-associated miR-21/29) is an
established node of muscle–heart cross-talk [4,5], motivating an explicitly
evidence-anchored regulatory analysis rather than database prediction alone.

Our overall design is summarised in Figure 1. We therefore (i) defined the shared
HF–muscle signature with statistical controls for the overlap itself; (ii) tested
whether it is separable from chronological ageing and HF etiology; (iii) mapped
upstream regulation using AGO-CLIP-supported miRNA evidence; and (iv) localised hub
genes at single-cell resolution. We emphasise at the outset that the muscle data
represent ageing muscle (old vs young) rather than clinically phenotyped
sarcopenia, and that the study is
hypothesis-generating.

---

## Results

### A statistically significant shared dysregulated gene set
Independent limma analyses (BH<0.05, |log₂FC|>0.378) yielded 921 HF and 610
muscle DEGs; their intersection (59 genes) significantly exceeded the random
expectation (hypergeometric p=1.0×10⁻³; 5,000-permutation p=8×10⁻⁴). Direction of
change among shared genes was **not** more concordant than chance (26/59
same-direction; binomial p=0.87), so we treated direction-concordance as a
*selection criterion* for a consistent panel rather than a discovered property.
Shared genes were enriched for extracellular matrix (GO-CC p_adj=1.3×10⁻⁴),
IL-4/IL-13 signalling (Reactome) and spanned ECM/fibrosis (LUM, SFRP1, MFAP4, DPT,
FBLN1, MYH11), NAD/mitochondrial (NAMPT, MTFP1, PPM1K) and mTOR (EIF4EBP1) genes.

### A parsimonious exploratory signature: strong in HF, modest in muscle
We used WGCNA disease modules (HF green r=−0.88; muscle turquoise r=−0.83) and six
ML selectors (≥3-vote concordance in both tissues) only to derive a *parsimonious
exploratory panel* of 11 hub genes, not as a primary discovery. These hubs were a
low-effect set: they were only modestly above a label-permutation null (p=0.075)
and were retained in the shared set at |log₂FC|>0.378 (11/11) but vanished at the
1.5-/2-fold cutoffs (0/11; threshold-sensitivity table, Supplementary Table S1),
i.e. the panel reflects many small-effect genes rather than a few strong markers.
A directional
signature score over the 8 hubs present across all cohorts (of the 11; 3 were
absent from the GSE5406 platform) discriminated HF in an independent
cohort (GSE5406 AUC 0.92, 95% CI 0.81–1.0). In ageing muscle, across **five
independent cohorts (N=134; GSE1428, GSE25941, GSE9103, GSE38718, GSE28392)** the
random-effects pooled AUC was **0.76 (95% CI 0.61–0.86; I²=0%)** — consistent and
above chance, with individual cohorts 0.64–0.83 (forest plot, Fig. 2C). Internal
and nested-CV estimates (0.96–1.0) are reported as apparent performance only.

### The HF component is robust to age and etiology
Although HF patients were older than controls (median 58 vs 52, p=1.4×10⁻⁵), age
alone discriminated weakly (AUC 0.64), whereas the signature retained AUC 0.95
after age adjustment and remained significant controlling for age (p=1.8×10⁻¹⁸).
A 4-gene non-ageing sub-signature (PITPNM1, MAGED2, MTFP1, FBLN1) preserved HF
discrimination (AUC 0.91) but lost muscle discrimination (0.54–0.70), indicating
the HF signal is robust to age whereas the muscle signal is more ageing-driven.
The signature performed equally in ischemic (AUC 0.96) and dilated (0.97)
cardiomyopathy, and the two etiologies were 92% concordant transcriptome-wide
(r=0.92), justifying their pooling for end-stage HFrEF.

### Upstream regulation converges on miR-15/16–CCND1
AGO-CLIP-supported miRNA→hub edges (ENCORI clipExpNum≥3; 257/436 also in
miRTarBase/TarBase) were dominated by CCND1 (72% of edges; a long-3′UTR target),
with the canonical miR-15/16/195/497 family — known regulators of cardiomyocyte
cell-cycle withdrawal [4,5] — the top multi-hub regulators (clipExpNum up to 67 on CCND1). Because clipExpNum partly tracks 3′UTR length and
study coverage, and the CLIP data derive from non-cardiac/muscle cell lines, we
present this network as hypothesis-generating. In a 3′UTR-length-matched null
(four matched control genes per hub), the hub panel as a whole was **not** enriched
for high-clipExp miRNAs (median 5.5 vs 4.0; paired Wilcoxon p=0.48), confirming
that aggregate edge counts largely reflect transcript length — with two
exceptions, CCND1 (204 vs 8.8) and NAMPT (51 vs 1), which substantially exceeded
their length-matched controls and which we therefore foreground as the
biologically interpretable nodes. ChEA3 added cardiac/fibrosis TFs (GATA4, SMAD1,
NKX2-5).

### Single-cell localisation
In 1.9M heart and 0.12M muscle cells (CELLxGENE Census [31]), NAMPT localised to
myeloid cells in both tissues (heart CD14⁺ monocytes 86%; muscle macrophages 67%)
and FBLN1 to stroma (cardiac fibroblasts 40%; tendon cells 94%), tying the shared
NAD/immune and ECM programmes to consistent cell lineages. Immune-deconvolution
(ssGSEA) showed prominent inflammatory remodelling in HF (17/22 cell types) and
modest changes in muscle (5/22).

---

## Discussion

We set out to test, rather than assume, whether failing myocardium and ageing
skeletal muscle share a definable transcriptional programme, and whether such an
overlap is separable from the chronological ageing common to both populations.
Three findings frame the interpretation.

First, the overlap is real but bounded. The 59-gene intersection significantly
exceeded chance (p=1.0×10⁻³), yet the shared genes were not preferentially
co-directional (concordance at chance) and were functionally generic —
extracellular-matrix/fibrosis, NAD/immune-metabolic and mTOR programmes that are
mobilised by tissue stress in many organs. We therefore frame the result as a
shared *stress-remodelling* programme rather than a bespoke organ-to-organ
messenger system, and we used direction-concordance as a panel-selection criterion
rather than presenting it as a discovered property.

Second, and most importantly for the ageing-confound concern, the HF arm of the
signature is robust to age: although HF patients were modestly older, age alone
discriminated weakly (AUC 0.64), whereas the signature retained AUC 0.95 after age
adjustment and a 4-gene non-ageing subset still reached 0.91. The muscle arm, by
contrast, weakened when ageing genes were removed, indicating that the shared
signal is driven by genuine HF biology on the cardiac side and by ageing biology
on the muscle side — a nuance that reframes the axis as "HF ↔ muscle ageing"
rather than HF ↔ a sarcopenia-specific lesion. Across five independent muscle
cohorts the signature reproducibly tracked muscle ageing (pooled AUC 0.76, I²=0%),
giving the muscle association a stable, if modest, footing.

Third, the regulatory and single-cell layers offer concrete, testable hypotheses.
Predicted upstream control converged on the canonical miR-15/16/195/497→CCND1
cell-cycle axis — a family with established roles in both cardiac remodelling and
muscle — although we caution that the CLIP evidence is dominated by CCND1 (a
long-3′UTR target heavily sampled in cancer CLIP) and is not tissue-specific.
Single-cell mapping localised the NAD/immune hub NAMPT to myeloid cells and the
ECM hub FBLN1 to stroma in both tissues, suggesting that the shared programme is
carried by analogous non-myocyte/non-myofibre compartments — a prediction that is
directly testable in cardiac and muscle single-cell or spatial datasets, and that
points to myeloid NAD metabolism and stromal ECM remodelling as candidate nodes of
a systemic cardio-muscle axis.

Clinically, the most provocative gap is that our HF data are end-stage explanted
HFrEF, whereas muscle wasting is most tightly linked to HFpEF in older adults. In
our GEO survey we did not identify a suitable bulk HFpEF *myocardial* case–control
series — the HFpEF datasets returned were blood, PBMC or epicardial-adipose —
consistent with HFpEF ventricle rarely being biopsied; this is a field-wide data
gap rather than a readily remediable omission. The natural next steps are
clinically phenotyped sarcopenic
muscle (e.g. the multi-omics sarcopenia cohort GSE226151, N=60) and any future
HFpEF myocardial series. We present this work as a hypothesis-generating, fully
reproducible in-silico map rather than a validated biomarker.

## Limitations
1. Muscle data are ageing muscle (old vs young), not clinically phenotyped
   sarcopenia; conclusions concern HF × muscle ageing.
2. HF cohort is end-stage explanted HFrEF (ischemic+dilated); HFpEF — most linked
   clinically to elderly muscle wasting — is absent.
3. Shared genes are functionally generic; direction-concordance is at chance and
   used only as a panel-selection criterion.
4. Muscle association, though reproducible across five cohorts (pooled AUC 0.76),
   is modest and individual cohorts are small (n=20–36); the muscle arm is partly
   ageing-driven (non-ageing subset AUC 0.54–0.70).
5. CLIP evidence is from non-target cell lines and CCND1-dominated; in a
   length-matched null the panel was not enriched beyond 3′UTR length (p=0.48,
   only CCND1/NAMPT exceeding control), so the network is hypothesis-generating,
   not context-specific proof.
6. Entirely in-silico; no functional validation.

## Methods
All parameters are fixed in `scripts/00_config.R`; the full pipeline (`scripts/00–17`,
`run_all.R`) is deterministic (seed 42).

**Datasets and groups.** HF discovery: GSE57338 left ventricle [9] (177 failing vs
136 non-failing; "heart failure:yes/no"). HF external validation: GSE5406 [8] (194
failing vs 16 donor). Ageing-muscle discovery: GSE8479 [6] (25 old vs 26 young; "Sample
Group" O vs Y, post-exercise samples excluded). Muscle validation (five cohorts):
GSE1428 [7] (description Older/Young), and GSE25941, GSE9103, GSE38718, GSE28392
[10] (Old/Young; baseline/sedentary samples for the exercise studies). Case =
failing/old, control = non-failing/young throughout.

**Preprocessing.** Series matrices were obtained from GEO [11] via GEOquery [12].
Expression was log₂-transformed when on a linear scale (GEO2R heuristic), probes mapped to symbols (Affymetrix "gene_assignment"
2nd token; Illumina GPL2700 via RefSeq→org.Hs.eg.db), collapsed gene-wise by the
max-mean probe, and quantile-normalised (`limma::normalizeBetweenArrays`). Cohorts
are kept separate (no cross-platform merging); each is z-scored within itself.

**DEG / WGCNA / shared set.** Per-tissue `limma` [13] (BH<0.05, |log₂FC|>0.378).
Signed WGCNA [14] on the top-5000-MAD genes (soft power by scale-free R²>0.85 → HF 12, muscle
14; minModuleSize 30, mergeCutHeight 0.25, deepSplit 2); disease module = strongest
module–trait correlation. Shared set = DEG_HF∩DEG_SAR (overlap significance by
hypergeometric test and 5,000-label permutations; direction-concordance by binomial
test vs the marginal up/down rate). Direction-concordant genes formed the panel
input.

**Hub panel and diagnostics.** Six selectors (LASSO [15], random forest [16],
SVM-RFE, Boruta [17], XGBoost [18], glmBoost [19]) were run in each discovery cohort; genes selected by ≥3
methods in *both* tissues defined the 11-hub exploratory panel (specificity vs a
200-label-permutation null reported). The diagnostic readout is a
platform-independent **directional signature score** = Σ (within-cohort z-scored
expression × discovery-derived sign); it contains no parameters fitted on
validation data. Discrimination = AUC with DeLong 95% CI (`pROC` [20]); HF additionally with
nomogram, 200×-bootstrap calibration and decision-curve analysis. Internal and
5-fold×20 nested-CV AUC (signs re-derived per training fold) are reported as
*apparent* performance only; external cohorts carry all inference.

**Muscle meta-analysis.** Per-cohort AUCs (DeLong CI) were logit-transformed and
pooled by random-effects meta-analysis (`metafor::rma` [21], REML; SE from the CI),
with I² for heterogeneity; the pooled estimate is back-transformed.

**Age / etiology robustness.** In GSE57338 (continuous age available) we compared
age-alone vs signature AUC, the age-residualised signature AUC, and the signature
coefficient in a logistic model adjusting for age; we recomputed AUC separately for
ischemic and dilated cardiomyopathy and correlated their per-gene log₂FC. A
non-ageing sub-signature was defined by removing hubs overlapping MSigDB
senescence/ageing sets.

**Regulatory network and bias handling.** miRNA→hub interactions from ENCORI/
starBase [25] were retained at clipExpNum≥3 (≥3 supporting AGO-CLIP experiments) and
cross-checked against miRTarBase [26]/TarBase [27]. Because clipExpNum co-varies
with 3′UTR length and study coverage and the CLIP libraries derive from
non-cardiac/muscle cell lines, per-hub edge counts, the edge–3′UTR-length
correlation and a 3′UTR-length-matched null were computed; the network is
interpreted as hypothesis-generating. Upstream TFs were obtained from ChEA3 [28]
(Integrated meanRank) and functional enrichment from clusterProfiler [22]
(GO/KEGG/Reactome) with STRING [24] PPI. A supplementary exploratory
drug-repurposing analysis used DSigDB/Enrichr [29,30].

**Single-cell and immune.** Hub expression was summarised across CZ CELLxGENE
Census v2023-12-15 (human heart 1.78M cells; musculature 0.12M) by streaming the
11-gene columns and accumulating per-cell-type mean log1p and percent-expressing
(cell types with ≥200 cells). Immune infiltration by ssGSEA (`GSVA` [23]) on canonical
immune-cell signatures.

## Data & code availability
All GEO accessions listed above. The complete deterministic pipeline
(`scripts/00..17` + `run_all.R`) and intermediate results are publicly available at
https://github.com/candicewu0515/HF-muscle-ageing-signature.

## References

**Clinical / biological background**
1. Anker SD, Ponikowski P, Varney S, et al. Wasting as an independent risk factor for mortality in chronic heart failure. *Lancet*. 1997;349(9058):1050–1053.
2. Tacke M, Ebner N, Boschmann M, Jarius A, Valentova M, Fülster S, et al. Resting energy expenditure and the effects of muscle wasting in patients with chronic heart failure (SICA-HF). *J Am Med Dir Assoc*. 2013;14(11):837–841. PMID 24094897; doi:10.1016/j.jamda.2013.08.008.
3. Cruz-Jentoft AJ, Bahat G, Bauer J, et al. Sarcopenia: revised European consensus on definition and diagnosis (EWGSOP2). *Age Ageing*. 2019;48(1):16–31.
4. Porrello ER, Mahmoud AI, Simpson E, et al. Regulation of neonatal and adult mammalian heart regeneration by the miR-15 family. *Proc Natl Acad Sci USA*. 2013;110(1):187–192. PMID 23248315; doi:10.1073/pnas.1208863110.
5. Porrello ER, Johnson BA, Aurora AB, et al. MiR-15 family regulates postnatal mitotic arrest of cardiomyocytes. *Circ Res*. 2011;109(6):670–679. PMID 21778430; doi:10.1161/CIRCRESAHA.111.248880.

**Datasets**
6. Melov S, Tarnopolsky MA, Beckman K, Felkey K, Hubbard A. Resistance exercise reverses aging in human skeletal muscle. *PLoS One*. 2007;2(5):e465. PMID 17520024; doi:10.1371/journal.pone.0000465. (GSE8479)
7. Giresi PG, Stevenson EJ, Theilhaber J, et al. Identification of a molecular signature of sarcopenia. *Physiol Genomics*. 2005;21(2):253–263. PMID 15687482; doi:10.1152/physiolgenomics.00249.2004. (GSE1428)
8. Hannenhalli S, Putt ME, Gilmore JM, et al. Transcriptional genomics associates FOX transcription factors with human heart failure. *Circulation*. 2006;114(12):1269–1276. PMID 16952980; doi:10.1161/CIRCULATIONAHA.106.632430. (GSE5406)
9. Liu Y, Morley M, Brandimarto J, et al. RNA-Seq identifies novel myocardial gene expression signatures of heart failure. *Genomics*. 2015;105(2):83–89. PMID 25528681; doi:10.1016/j.ygeno.2014.12.002. (GSE57338)
10. Human skeletal-muscle ageing series GSE25941, GSE9103, GSE38718, GSE28392 (NCBI GEO; cited by accession).

**Methods, tools and databases**
11. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets—update. *Nucleic Acids Res*. 2013;41(D1):D991–D995.
12. Davis S, Meltzer PS. GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics*. 2007;23(14):1846–1847.
13. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res*. 2015;43(7):e47.
14. Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis. *BMC Bioinformatics*. 2008;9:559.
15. Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent. *J Stat Softw*. 2010;33(1):1–22. (glmnet)
16. Breiman L. Random forests. *Mach Learn*. 2001;45(1):5–32.
17. Kursa MB, Rudnicki WR. Feature selection with the Boruta package. *J Stat Softw*. 2010;36(11):1–13.
18. Chen T, Guestrin C. XGBoost: a scalable tree boosting system. *Proc 22nd ACM SIGKDD*. 2016:785–794.
19. Hothorn T, Bühlmann P, Kneib T, et al. Model-based boosting 2.0. *J Mach Learn Res*. 2010;11:2109–2113. (glmboost)
20. Robin X, Turck N, Hainard A, et al. pROC: an open-source package for R and S+ to analyze and compare ROC curves. *BMC Bioinformatics*. 2011;12:77.
21. Viechtbauer W. Conducting meta-analyses in R with the metafor package. *J Stat Softw*. 2010;36(3):1–48.
22. Wu T, Hu E, Xu S, et al. clusterProfiler 4.0: a universal enrichment tool for interpreting omics data. *Innovation (Camb)*. 2021;2(3):100141.
23. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7.
24. Szklarczyk D, Kirsch R, Koutrouli M, et al. The STRING database in 2023. *Nucleic Acids Res*. 2023;51(D1):D638–D646.
25. Li J-H, Liu S, Zhou H, et al. starBase v2.0: decoding miRNA–ceRNA, miRNA–ncRNA and protein–RNA interaction networks from large-scale CLIP-Seq data. *Nucleic Acids Res*. 2014;42(D1):D92–D97. (ENCORI)
26. Huang H-Y, Lin Y-C-D, Cui S, et al. miRTarBase update 2022. *Nucleic Acids Res*. 2022;50(D1):D222–D230.
27. Karagkouni D, Paraskevopoulou MD, Chatzopoulos S, et al. DIANA-TarBase v8. *Nucleic Acids Res*. 2018;46(D1):D239–D245.
28. Keenan AB, Torre D, Lachmann A, et al. ChEA3: transcription factor enrichment analysis by orchestration of multiple libraries. *Nucleic Acids Res*. 2019;47(W1):W212–W224.
29. Kuleshov MV, Jones MR, Rouillard AD, et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res*. 2016;44(W1):W90–W97.
30. Yoo M, Shin J, Kim J, et al. DSigDB: drug signatures database for gene set analysis. *Bioinformatics*. 2015;31(18):3069–3071.
31. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. *bioRxiv*. 2023; doi:10.1101/2023.10.30.563174. (CELLxGENE Census)
