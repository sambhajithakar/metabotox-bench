# MetaboTox-Bench — mechanism gold set, v0.1

**Companion to:** [metabolite-toxicity-ooc-gap-analysis.md](metabolite-toxicity-ooc-gap-analysis.md)
**Data file:** [metabotox_bench_gold_set_v0.2.csv](metabotox_bench_gold_set_v0.2.csv) — **use v0.2, not v0.1**
**References:** [gold_set_references_draft.csv](gold_set_references_draft.csv)
**Built:** 8 October 2026 · curated to v0.2 overnight 8–9 October · 46 rows · 12 columns
**Target venue:** *Journal of Chemical Information and Modeling* (ACS, hybrid — subscription route, no APC). Backup: *Archives of Toxicology*.

> **v0.2 changed six rows and removed one control.** Reference collection turned up evidence that contradicted the original coding in several places. v0.1 is kept only for provenance; do not model on it. See §3b.

**Locked project decisions**

- **Final title:** *Metabolite-Resolved Multi-Task Deep Learning for Cross-Organ Off-Target Toxicity: A Multiple-Instance Framework with Calibrated Uncertainty and a Public Benchmark*
- **Journal:** *Briefings in Bioinformatics*. Framing is therefore **methods-first**: the MIL formulation and this benchmark lead the paper; toxicology is the application domain. Confirm the APC policy on OUP's own page on submission day, and check for an institutional read-and-publish agreement.
- Organ-on-a-chip is **validation and interpretation**, not a training substrate, and stays out of the title.

---

## 1. What this set is, and what it is not

This is **not** the training set. It is the **mechanism evaluation set** — the part of MetaboTox-Bench that tests the paper's actual claim.

The claim is not "we predict toxicity better." It is: *given a parent drug, the model names which metabolite carries the hazard and which organ is injured.* Standard benchmarks cannot test that, because their labels are indexed by parent drug. These 46 rows can, because for each one the field already knows the answer.

**Training labels still come from the public corpora** listed in §5 of the gap analysis (Tox21, TDC, DILIrank/LTKB, SIDER, DILImap). This set is held out entirely and never trained on.

## 2. Inclusion criteria

A drug enters only if the toxic species, the bioactivating step, and the affected organ are each established in the toxicology literature well enough to act as ground truth. That is why the set is small; it is meant to be small and correct, not large.

## 3. Composition

| Case type | n | Role in evaluation |
|---|---|---|
| `positive_metabolite` | 39 | The metabolite is the hazard. Attribution must name it. |
| `inverse` | 1 | **Parent is the hazard; the metabolite is safer.** Terfenadine (→ fexofenadine). **Only one left after curation — this is the set's main weakness.** |
| `pgx_inverse` | 1 | **Failure to metabolise causes toxicity** — perhexiline in CYP2D6 poor metabolisers. |
| `negative_metabolite` | 2 | Real toxicity, **not** metabolite-driven — bosentan (parent BSEP inhibition) and tolcapone (parent uncoupling). |
| `negative_control` | 1 | Near-identical scaffold to a toxic drug, no toxicity — entacapone. |
| `parent_and_metabolite` | 1 | Both parent and metabolite hazardous — astemizole. |
| `pgx` | 1 | Toxicity set by metabolite accumulation — mercaptopurine/azathioprine (TPMT, NUDT15). |

**28 of 46 are cross-organ** (`cross_organ = Y`): the metabolite is formed in one organ and injures another. Those rows are the ones that actually test the paper's central claim. 30 rows are `high` confidence, 16 `medium` after curation.

## 3b. What curation changed, and why v0.1 must not be used

Collecting a primary reference for every row surfaced evidence that contradicted the original coding. Six rows were corrected and one control was lost. This is recorded here and in the manuscript rather than quietly fixed, because the corrections change what the set tests.

| Row | Was | Now | Evidence |
|---|---|---|---|
| **Astemizole** | `inverse` — metabolism detoxifies | `parent_and_metabolite` — both hazardous | Desmethylastemizole blocks hERG with potency comparable to the parent (Zhou et al. 1999, doi:10.1111/j.1540-8167.1999.tb00264.x) and torsade has been attributed to the metabolite specifically (PMID 8917271). **No longer usable as an inverse control.** |
| **Methoxyflurane** | fluoride + oxalic acid | fluoride **plus dichloroacetic acid** | Kharasch et al. 2006, *Anesthesiology* 105(4):726–736 — fluoride alone caused diuresis but **not** significant tubular necrosis; necrosis required fluoride together with co-formed DCAA. **Oxalic acid is not implicated.** |
| **Tolcapone** | `positive_metabolite` | `negative_metabolite` | Established mechanism is mitochondrial uncoupling by the **parent** nitrocatechol; no study found establishing a nitroreduction metabolite as mediator. The entacapone pair now tests **parent-level** discrimination. |
| **Zomepirac** | liver, anaphylaxis | **kidney** added | The strongest causal evidence for acyl-glucuronide mediation (Iwamura et al. 2016) demonstrates acute kidney injury, not hepatotoxicity. Confidence lowered. |
| **Primaquine** | hydroxylamine | **5,6-orthoquinone**, formation partly in erythrocyte | Current mechanistic work implicates the orthoquinone downstream of CYP2D6-derived 5-hydroxyprimaquine; activation occurs partly in the red cell. **Strengthens the tissue-conditioning test.** |
| **Capecitabine** | liver, tumour | **skin** added as forming organ | Thymidine phosphorylase converts 5′-DFUR to fluorouracil locally in skin in hand-foot syndrome. Formation and injury co-localise — **a point in favour of the organ-resolved design.** |
| **Cyclophosphamide** | acrolein, high confidence | acrolein, **medium** | Acrolein-only attribution contested by an acrolein-independent urothelial apoptotic wave (PMID 24353901) and a chloroacetaldehyde contribution (PMID 31598736). The mesna-rescue argument is weaker than it looked. |

Confidence was also lowered on **nefazodone** (no primary study for the quinone imine as mediator, reviews only), **leflunomide** (hepatic P450s *attenuate* cytotoxicity, so teriflunomide-as-toxicant is not clean) and **carbamazepine** (HLA alleles dominate the endpoint; metabolite causation is hypothesis). **Sudan I** was refined to the benzenediazonium ion, with a note that Sudan I is also a CYP1A1 *inducer* — a different role from being a substrate, which the current schema conflates and should separate.

**The open weakness:** only one clean inverse control (terfenadine) survives. One row cannot establish that a model has learned the protective-metabolism case. Expanding this class is the highest-value addition to v0.3; candidate cases to verify include cisapride (norcisapride far less active) and other drugs whose marketed metabolite is a separate, safer product. Do not add any of them without a primary reference.

---

## 4. Why the non-positive cases matter more than their count suggests

A model that has learned "reactive functional group present → blame the metabolite" will score well on the 40 positives and **fail every one of the other six**. That is exactly the failure mode a reviewer will look for, so these are the rows that decide whether the paper holds up:

- **Terfenadine and astemizole** — metabolism *removes* the hazard. A model that blames fexofenadine is wrong in the most visible possible way.
- **Perhexiline** — toxicity from *not* metabolising. Breaks the assumption that metabolism creates hazard.
- **Bosentan** — real hepatotoxicity that must **not** be attributed to any metabolite.
- **Entacapone vs tolcapone** — a matched pair on near-identical scaffolds, one hepatotoxic and one not. Structure-only models tend to flag both. This is the single sharpest discriminating pair in the set, and it is the same pair used as the real-world reference in TransTox.
- **Amiodarone** — a toxic metabolite (desethylamiodarone) that is **not electrophilic**; it accumulates and causes phospholipidosis. Guards against a model that has only learned reactive-group alerts.

Three further rows are deliberately awkward and worth keeping for that reason: **clozapine** and **nevirapine** (bioactivation partly *extrahepatic* — neutrophil MPO and skin sulfotransferase, so "liver forms it" is wrong); **benzene** (two-step, two-organ: hepatic phenol then marrow MPO oxidation); **primaquine** and **dapsone** (target organ is the erythrocyte, which no standard organ panel includes — a prompt to settle the organ vocabulary before modelling).

## 5. Planned evaluation protocol — fix this before any modelling

Per §7 step 3 of the gap analysis, the protocol is frozen in advance so that a flat AUROC cannot be rationalised after the fact.

**Primary endpoints**
1. **Top-1 metabolite attribution accuracy** on the 40 positives — does the highest-attributed instance match the known toxic species?
2. **Organ-pair accuracy** on the 28 cross-organ rows — are both forming and affected organ correct?
3. **Control compliance** — the 6 inverse/negative rows must *not* be attributed to a metabolite. Report this as a hard pass/fail per row, not an average.

**Secondary**
4. Calibration (ECE, Brier) on parent-level risk, with the applicability-domain flag reported.
5. AUROC/AUPRC on standard splits, versus parent-only multi-task baselines (FATE-Tox, MT-Tox) and the rule-based chaining baseline — reported as **comparability, not superiority**.

**Splits:** scaffold split *and* a metabolite-disjoint split, so that no metabolite structure appears in both train and test.

## 6. Known limitations — state these in the paper

1. **References exist but 14 are incomplete.** All 46 rows now carry a verified reference (37 primary mechanistic studies, 8 reviews, 1 CPIC guideline) in `gold_set_references_draft.csv`. For 32 rows full bibliographic metadata was resolved; **14 carry a verified title and indexer URL only** and must be completed against a reference manager before release. Four year ambiguities (troglitazone, amodiaquine, CPIC, bosentan — online-first versus print issue) are flagged in that file and were deliberately left as the index gave them rather than corrected by guess.
2. **Selection bias toward the well-understood.** These are the cases the field solved, which are not representative of the idiosyncratic reactions that matter most clinically. The set measures mechanistic attribution, not clinical coverage — say so explicitly.
3. **Withdrawn-drug enrichment.** Many rows are withdrawn or black-boxed drugs, so they are over-represented in historical training corpora. Check train/test leakage against DILIrank and LTKB before reporting.
4. **Organ vocabulary is unsettled.** Erythrocyte, bone marrow, urothelium and "systemic autoimmune" do not fit a standard five-organ panel. Decide the ontology (and how to handle `systemic`) before building label tensors.
5. **13 rows are `medium` confidence.** Either upgrade them with references or report results with and without them.
6. **Hydralazine and procainamide are systemic, not organ-localised.** Either define a `systemic` class or drop both.

## 7. Next steps

1. Attach a primary reference to each of the 46 rows; drop or demote anything that will not support one.
2. Resolve the open question from §7b of the gap analysis — read the ACS Ames-mutagenicity MIL paper and confirm what its instances are, since the novelty wording depends on it.
3. Fix the organ ontology (§6.4), then freeze the label schema.
4. Run the metabolite generators (BioTransformer, SyGMA, MetaPredictor) over all 46 parents and measure **recall of the known toxic species** — if a generator never proposes NAPQI, no downstream model can attribute it. This is a hard upstream ceiling and should be measured and reported early.
5. Assemble the training corpus and build the two splits.

---

## v0.3 — 12 October 2026

**Use `metabotox_bench_gold_set_v0.3.csv`.** v0.2 and v0.1 are kept for provenance only.

One change from v0.2, found by the structure-provenance protocol in `../run1/`
(see `run1_README.md` §4), not by re-reading the references:

**Terbinafine.** The toxic species was recorded as *TBF-A (7,7-dimethylhept-2-ene-4-ynal)*.
That name is wrong. Converting it with OPSIN, a deterministic IUPAC name parser, produced a
structure with an **isopropyl** terminus, whereas TBF-A retains terbinafine's *tert*-butyl
group. The correct identity is **(E)-6,6-dimethylhept-2-en-4-ynal**.

Source: Barnette DA, Davis MA, Flynn N, Pidugu AS, Swamidass SJ, Miller GP. Comprehensive
Kinetic and Modeling Analyses Revealed CYP2C9 and 3A4 Determine Terbinafine Metabolic
Clearance and Bioactivation. *Biochem Pharmacol* 2019;170:113661. PMID 31605674;
PMCID PMC6905088; doi:10.1016/j.bcp.2019.113661.

Also updated from the same reference: `bioactivating_enzyme` is now
`CYP2C9;CYP3A4 (N-dealkylation)` rather than the unspecific `CYP_N_dealkylation`, and
`citation_status` is `verified_primary`.

**Why this is worth recording as a method note.** Running a curated *name* through a
deterministic converter is a cheap integrity check on a mechanism set: a wrong name becomes
a visibly wrong structure instead of being absorbed silently. It also mattered for the
result — with the wrong structure, terbinafine would have been scored a permanent
generator miss instead of the one species that search depth recovers.

**Structure coverage.** Verified parent-and-species structure pairs exist for 18 of the 46
entries; see `../run1/gold_set_structures_v0.1.csv` for the table and its per-row
provenance, and `run1_README.md` §3 for the two admissible routes. The 28 entries without
a pair are mostly reactive intermediates with no reachable public structure record, and
that non-retrievability is reported as a finding rather than treated as a smaller
denominator.
