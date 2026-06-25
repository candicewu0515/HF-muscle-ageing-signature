# Peer review + added robustness evidence — honest verdict

Self-audit (`scripts/self_check.R`) + adversarial peer review + added analyses
(`scripts/12_robustness.R`, `13_nestedcv.R`). Outputs in `results/13_robust/`.

## What the added evidence shows (the uncomfortable truth)

| Test | Result | Verdict |
|------|--------|---------|
| Overlap of HF & muscle DEGs | 59 observed vs ~28 expected, **hypergeometric p=1.0e-3** | ✅ real — shared dysregulation is significant |
| Direction concordance of overlap | 26/59 same-direction, **binomial p=0.87** (null rate 0.51; observed 44% < chance) | ❌ **at chance** — "cross-tissue direction consistency" is NOT a discovered signal |
| ML hub set vs label-permutation null | observed 9 shared-selected, null mean 6.5, **p=0.075** | ⚠️ only borderline above chance |
| DEG-threshold sensitivity | hubs retained 11/11 at logFC 0.378, **0/11 at 0.585 / 1.0** | ❌ entire hub set is low-fold-change & threshold-dependent |
| Aging/senescence overlap (MSigDB) | **16/26 core genes** in aging sets, **p=1.8e-11**; 7/11 hubs | ⚠️ signature is **majority aging biology** |
| CLIP network composition | **CCND1 = 72% of all 436 edges**; rho(edges,3'UTR len)=0.54 | ❌ single-gene-dominated; clipExpNum tracks 3'UTR length / study popularity, not tissue-specific binding |
| Nested-CV AUC (signs re-derived per fold) | HF 0.96 / SAR 1.0 — but genes pre-selected on full cohort | internal numbers are **apparent (optimistic)**; external cohorts are the only leakage-free estimate |

## Five attacks that survive (must be addressed, not hidden)

**A. Circularity.** Internal AUC≈1.0 is the fingerprint of selection on the same
data. Re-deriving signs per fold barely moves it (0.96–1.0) because the *genes*
were chosen on this cohort. → Lead with **external** AUC only (HF 0.92, SAR
0.64–0.72); label all internal/nested numbers "apparent".

**B. Aging ≠ sarcopenia.** GSE8479 is old-vs-young healthy muscle, no clinical
sarcopenia phenotype. 62% of the core overlaps senescence/aging sets (p=1.8e-11);
HF patients skew old. The shared signal is plausibly **chronological aging**, not
a sarcopenia-specific mechanism. → Re-label throughout as "HF × **aging skeletal
muscle / muscle wasting**", state age-confound as central limitation.

**C. Missing tests now supplied.** Overlap significant (p=1e-3); concordance NOT
(p=0.87); ML hubs borderline (p=0.075). Report all three honestly.

**D. CLIP artifact.** 72% of edges are CCND1 (long-3'UTR oncogene, saturated in
pan-cancer AGO-CLIP; the CLIP data is from cancer cell lines, not heart/muscle).
→ Demote "CLIP-anchored hard evidence" to "hypothesis-generating"; foreground the
single well-validated **miR-15/16/195/497 → CCND1** axis honestly; drop the
"N AGO-CLIP experiments = context-specific proof" framing.

**E. Small n + threshold.** SAR external n=22/36, AUC CI crosses 0.5; hubs vanish
at 1.5-fold. → Report SAR as exploratory; add threshold-sensitivity table (done);
ideally pool more muscle cohorts for a meta-AUC.

## Honest revised story (what actually holds)

> Heart failure and **aging skeletal muscle** share a significantly-overlapping,
> **senescence-enriched** dysregulated gene set. A modest-effect 8-gene subset
> behaves as a diagnostic signature — robust in HF (external AUC 0.92), weaker in
> aged muscle (0.64–0.72). The shared genes converge on ECM/fibrosis, NAD–immune
> metabolism (NAMPT, localized to myeloid cells in both tissues) and mTOR
> (EIF4EBP1); their predicted upstream regulation is dominated by the canonical
> **miR-15/16 → CCND1** cell-cycle axis. The work is a hypothesis-generating
> in-silico map of a shared aging-associated cardio-muscle axis, not a validated
> sarcopenia biomarker.

## Added hard evidence (scripts/14) — partially rescues the critiques

- **Age confound (HF):** HF older than controls (median 58 vs 52, p=1.4e-5), but
  age alone gives AUC 0.643, while the signature gives 0.965 and **0.953 after
  age adjustment** (signature p=1.8e-18 controlling for age). → HF discrimination
  is **age-independent**, not chronological aging.
- **Non-aging sub-signature:** dropping the 7 aging hubs, the 4 non-aging hubs
  (PITPNM1, MAGED2, MTFP1, FBLN1) still give HF AUC **0.914** (vs 0.917 full).
  Muscle drops (0.54–0.70) → muscle signal is more aging-dependent (honest nuance).
- **HF subtype:** signature works in ischemic (0.961) and dilated (0.969) CM; the
  two etiologies are 92% concordant transcriptome-wide (rho=0.92) → pooling
  ICM+DCM is justified for end-stage HFrEF. (HFpEF still absent — limitation.)

## Verdict
Peer-review consensus: **as currently framed (sarcopenia-specific, direction-
consistent, CLIP-proven) → major revision / likely reject.** Reframed honestly
(aging-muscle, hypothesis-generating, external-only performance) it is a
defensible 0–3 journal submission. The added robustness analyses ARE the fix —
presenting them pre-empts every attack. Highest-value remaining work: (1) pool
more muscle cohorts for SAR meta-AUC; (2) length/abundance-matched CLIP null;
(3) age-adjusted DEG to separate aging from HF-specific signal.
