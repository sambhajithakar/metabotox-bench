# MetaboTox-Bench

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23402205.svg)](https://doi.org/10.5281/zenodo.23402205)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Data and code for:

> **Machine learning cannot exceed its inputs: rule-based metabolite enumeration bounds
> AI-driven organ toxicity attribution.**
> S. B. Thakar and P. N. Ghorpade, Advanced Bioinformatics, Data Science and Generative AI Unit,
> Aditya BioNova Analytics, Thakar Mala Road, At/Post Haroli, Taluka Shirol,
> District Kolhapur 416102, Maharashtra, India.
> Submitted to Toxicological Sciences, 2026.

Every number in the paper comes from a file in this repository. Nothing is estimated.

## What the paper measures

Pipelines that attribute organ toxicity to a metabolite all begin by *enumerating* candidate
metabolites. A species a generator never proposes cannot be blamed by any model downstream,
so enumeration coverage is a hard ceiling on the whole pipeline — and it is almost never
reported.

On 18 entries whose toxic species is established in the primary literature and whose
structure is publicly retrievable, under one matching routine and one denominator:

| Generator | Rule applications | Coverage | Mean candidate set |
|---|:--:|:--:|--:|
| SyGMa 1.1.0 | 2 | 4/18 = 22.2% | 18 |
| SyGMa 1.1.0 (documented default) | 3 | 7/18 = 38.9% | 82 |
| SyGMa 1.1.0 | 4 | 8/18 = 44.4% | 263 |
| GLORYx rule set, exhaustive | 1 | 7/18 = 38.9% | 14 |
| GLORYx rule set, exhaustive | 2 | **13/18 = 72.2%** | 170 |
| GLORYx rule set, exhaustive | 3 | **15/18 = 83.3%** | 1992 |

Full-InChIKey and skeleton matching agree at every setting, so every miss is a structural
absence rather than a stereochemical near miss. All four parent-driven control entries are
recovered throughout, which places every failure in the generative step rather than in the
matching.

Because both generators run on the same 18 entries, the comparison is paired. For the
GLORYx rule set at two applications against SyGMa at its documented default, 7 entries are
recovered by both, 5 by neither, 6 by the GLORYx rules alone and **none by SyGMa alone**
(exact McNemar on the 6 discordant pairs, p = 0.031). On these entries one rule set's
coverage is a superset of the other's, which is a stronger statement than the difference in
the two percentages.

**Two results.** At matched candidate-set cost, the choice of rule set is worth roughly
thirty points of attainable coverage. And three entries resist both rule sets at every depth
enumerated — methoxyflurane → dichloroacetic acid, leflunomide → teriflunomide, and
benzene → benzene oxide. Benzene is the diagnostic case: phenol, hydroquinone and
1,4-benzoquinone are all covered, but the arene oxide between benzene and phenol is not.
Rule sets encode net transformations and step over transient intermediates, and transient
intermediates are disproportionately the species that do the damage.

## Two limitations that bound every GLORYx number here

1. **GLORYx itself was not run.** GLORYx is FAME 3 site-of-metabolism prediction plus
   reaction rules, and the FAME 3 models are not publicly distributed. What is measured here
   is the **expressive ceiling of the published rule set**: all 260 rules applied
   exhaustively, no site filtering, no scoring. Since FAME 3 can only select among rule
   products, this is a generous **upper bound** on the tool, not a measurement of it.
2. **The two generators are not independent.** 178 of GLORYx's 260 rules derive from SyGMa,
   so the three shared failures are weaker evidence than agreement between unrelated tools
   would be.

## Layout

| Path | What it is |
|---|---|
| `metabotox_bench_gold_set_v0.3.1.csv` | **The current reference set** — 46 entries, 13 columns |
| `metabotox_bench_gold_set_v0.3.csv` | The version archived under the Zenodo DOI, kept for provenance |
| `metabotox_bench_gold_set_README.md` | Its documentation and change log |
| `gold_set_references_draft.csv` | One verified primary reference per entry |
| `run1/gold_set_structures_v0.1.csv` | Verified structure table, with provenance per row |
| `run1/build_structures.py` | Builds it, both provenance routes |
| `run1/run1_generator_recall.py` | SyGMa coverage harness |
| `run2/run2_gloryx_rules_recall.py` | GLORYx rule-set coverage harness, with its tripwire |
| `run2/gloryx_reactionrules.csv` | The 260 published GLORYx rules, verbatim — **GPL-3.0, see `run2/GLORYX_NOTICE.md`** |
| `run2/run2_results_d*.csv` | Per-entry results, 1 / 2 / 3 applications |
| `run2/run2_residual_d4.csv`, `_d5.csv` | Depth-4 and depth-5 test on the three residual failures |
| `run2/*.log` | Console logs, including every self-check and cap report |
| `run2/run2_summary.md` | Full run record: what was and was not run, both harness traps |
| `figures/make_figures.py` | Regenerates all three figures from the result files |

## Reproducing

```
python -m venv venv && . venv/bin/activate
pip install "setuptools<81"
pip install --no-build-isolation sygma       # see run1/run1_README.md for why
pip install rdkit py2opsin matplotlib
python run2/run2_gloryx_rules_recall.py      # see the script for depth arguments
python figures/make_figures.py
```

`figures/make_figures.py` reads every number from the result files rather than having them
typed in, and asserts that all four parent-driven control rows are recovered in every file it
reads, so a harness regression cannot be plotted silently.

## Two traps worth knowing about if you build on this

1. **GLORYx's SMIRKS need explicit hydrogens.** They are written for CDK and name hydrogen as
   a real atom, so on an implicit-hydrogen RDKit molecule almost every rule fails to match
   *without raising anything*. Our first run reported 38.9% at two applications instead of
   72.2%, with a mean candidate set of 4 instead of 170 — and the wrong number was plausible
   because it happened to equal SyGMa's default. `run2_gloryx_rules_recall.py` now refuses to
   print any coverage figure unless paracetamol yields at least ten products in one step with
   NAPQI among them.
2. **Keep the parent in the candidate bag.** SyGMa's `to_list()` includes the root. Removing
   the parent silently destroys the four parent-driven control rows while leaving the eighteen
   real rows untouched, so the control stops controlling and nothing looks wrong.

## Which version of the reference set to use

Use **`metabotox_bench_gold_set_v0.3.1.csv`**. It differs from the archived `v0.3` in two
ways and in no others — every mechanism, structure, organ and endpoint field is byte-identical
across all 46 rows:

- a new `primary_reference` column carrying the citation key of the primary reference for each
  entry, matching Table 5 of the paper;
- `citation_status` now reads `verified_primary` for all 46 rows, where `v0.3` still read
  `needs_primary_ref` for 45 of them.

`v0.3` was written before the reference curation was finished, so its `citation_status`
column contradicted the paper's Table 5. The paper states this discrepancy rather than
quietly correcting it, and `v0.3` is kept here unchanged so the archived record and this
repository can be compared.

## Relation to the Zenodo deposit

The archival record is **https://doi.org/10.5281/zenodo.23402205** (17 October 2026, 37
files). That DOI is what the paper cites and what should be cited in any reuse. This
repository is the working copy: it carries the same code and results plus `v0.3.1` of the
reference set and the 400 dpi TIFF figures, and it is where corrections land first. Where the
two disagree, the DOI is the record of what was submitted and this repository is the current
state.

## Licence

Our code and data: MIT (see `LICENSE`).
`run2/gloryx_reactionrules.csv` is third-party content under **GPL-3.0** — see
`run2/GLORYX_NOTICE.md` and `run2/GLORYX_LICENSE_GPL-3.0.txt`.
