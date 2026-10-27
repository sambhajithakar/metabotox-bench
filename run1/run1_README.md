# Run 1 — generator recall ceiling

**Built and executed 12 October 2026.** This is the experiment the submission checklist
put first, because it could redirect the whole project. It did.

**Why it was first:** if a metabolite enumerator never proposes the established toxic
species, no downstream model can attribute toxicity to it. Generator recall is a hard
ceiling on every attribution number in the paper. Measuring it costs almost nothing.

---

## 1. Result

SyGMa 1.1.0, phase-I depth 2 / phase-II depth 1 (the tool's documented default), on the
22-row verified structure table:

```
Rows in table:        22
Evaluated:            22
Parent-driven rows (recall not applicable): 4 -> Terfenadine, Bosentan, Tolcapone, Perhexiline
Recall denominator:   18

  sygma:
    recall (full InChIKey):       7/18 = 38.9%
    recall (skeleton, no stereo): 7/18 = 38.9%
    mean bag size:                82
    MISSED: Valproic acid, Cyclophosphamide, Terbinafine, Felbamate, Halothane,
            Methanol, Methoxyflurane, Leflunomide, Benzene, Benzene, Benzene
```

**Depth sweep.** Recall plateaus; the candidate set does not.

| Phase-I depth | Recall | Mean bag size | Recovered vs. depth 2 |
|---|---|---|---|
| 2 (default) | 7/18 = 38.9% | 82 | — |
| 3 | 8/18 = 44.4% | 263 | terbinafine → TBF-A |
| 4 | 8/18 = 44.4% | 657 | none |

Eightfold more candidates for one extra species. Depth is neither the explanation nor
the remedy.

**Two checks that make the number trustworthy:**

- **Full and skeleton recall are identical at every depth.** Not one miss is recovered by
  ignoring stereochemistry. The failures are absences, not near misses.
- **All 4 parent-driven controls are recovered.** The matching procedure finds a species
  when it is there, so the failures are generative, not an artefact of matching or of our
  structures.

**The failures are chemically coherent** — five classes the rule set structurally cannot
express:

| Class | Cases |
|---|---|
| Scaffold fragmentation | cyclophosphamide → acrolein; felbamate → atropaldehyde |
| Desaturation | valproic acid → 4-ene-VPA |
| Oxidative dehalogenation | halothane → trifluoroacetyl chloride; methoxyflurane → dichloroacetic acid |
| Heterocyclic ring-opening | leflunomide → teriflunomide |
| Small, non-drug-like substrates | methanol → formate; benzene → oxide, hydroquinone, 1,4-benzoquinone |

Benzene is the starkest: **SyGMa returns a single candidate for it.** Methanol, three.

## 2. What this did to the project

The checklist pre-committed to a decision at this gate: if recall came back weak, split
the paper and publish the benchmark with the recall diagnosis as the central finding,
rather than stacking a model on a broken upstream step. 38.9% is weak, so that is what
happened. See `../manuscript/00_READ_ME_FIRST.md` and `../manuscript/MANUSCRIPT_READY.md`.

## 3. How the structure blocker was broken

The gold set records toxic species by *name*; the harness needs SMILES. PubChem's REST
API, ChEMBL's API, Europe PMC and direct `curl` egress are all refused by this container's
proxy, and writing 92 SMILES from memory was not acceptable — a mis-drawn reactive
intermediate yields a recall number wrong in a direction nobody can audit.

Two routes turned out to work, and only these two were admitted:

**Route A — Wikipedia infobox with a recomputation check.** The page gives SMILES *and*
InChIKey. Recompute the InChIKey from the SMILES with RDKit and compare: agreement rules
out transcription error in either field, since an error in one would not reproduce the
other. **Twenty parent structures were admitted this way and all twenty matched exactly**,
stereochemical block included. Note that drugboxes render these fields to a fetcher but
chemboxes (benzene, hydroquinone, hydrazine, oxalic acid, dichloroacetic acid) do not —
which is why Route B is needed.

**Route B — OPSIN on a name taken from the gold set.** `py2opsin` is a deterministic IUPAC
name-to-structure converter, so the only judgement is the name, and the name comes from
the curated gold set rather than from recall. Nineteen species were admitted this way.
Validated end to end against Route A: *N*-(4-hydroxyphenyl)acetamide → OPSIN → RDKit
InChIKey reproduced Wikipedia's paracetamol key `RZVAJINKPMORJF-UHFFFAOYSA-N` exactly.

Anything outside these two routes was left out and reported as an exclusion. That leaves
**18 of 46 gold-set entries (39%) with a verified parent-and-species pair.** The other 28
are a finding in their own right: several canonical toxic species — the aristolactam
nitrenium ion, the amodiaquine quinone imine, the thiophene *S*-oxides — have no public
structure record reachable from here at all.

## 4. A gold-set error the protocol caught

The gold set recorded terbinafine's reactive aldehyde TBF-A as
*7,7-dimethylhept-2-ene-4-ynal*. OPSIN converted that faithfully, to a structure with an
**isopropyl** terminus — which is wrong, because TBF-A retains terbinafine's *tert*-butyl
group. The correct name is **(E)-6,6-dimethylhept-2-en-4-ynal** (Barnette DA, Davis MA,
Flynn N, Pidugu AS, Swamidass SJ, Miller GP. *Biochem Pharmacol* 2019;170:113661;
PMID 31605674; PMCID PMC6905088). Corrected in gold set **v0.3**, along with the
bioactivating enzymes, which that paper gives as CYP2C9 and CYP3A4.

This is the protocol working as designed: a deterministic converter turns a naming error
into a visible structural anomaly instead of absorbing it. It also matters for the result
— with the wrong name, terbinafine would have been scored a permanent miss instead of a
depth-3 hit.

## 5. Environment

The system Python here is 3.13, for which no RDKit wheel exists. That is why the obvious
install fails. What works:

```
python3.12 -m venv rdenv
./rdenv/bin/pip install rdkit "setuptools<81"
./rdenv/bin/pip install --no-build-isolation sygma
./rdenv/bin/pip install py2opsin          # needs java, which is present
```

Both extra flags are required: SyGMa's `setup.py` checks for RDKit InChI support that the
isolated build environment lacks, and it imports `pkg_resources`, removed in setuptools 81+.

**One real bug, which anyone reproducing this will hit:** SyGMa's `to_list()` sorts on a
score field that is `None` until `calc_scores()` has been called, and raises
`TypeError: '<' not supported between instances of 'NoneType' and 'NoneType'`. Call
`tree.calc_scores()` first. The metabolite score key is `SyGMa_score`.

## 6. Still to wire up — and this bounds the claim

**Only SyGMa was run.** This is the single biggest limitation of the result and it is
stated as such in the manuscript (§7.1): the honest claim is about *one rule-based
generator*, not about the field's ceiling.

- **BioTransformer** — not implemented. Needs the BioTransformer 3 jar, absent here;
  `java` *is* present, so this is a download-and-invoke job, not a porting job.
- **MetaPredictor** — not implemented. Needs model weights not available here.
- Neither stub returns a number; both raise `NotImplementedError`, so a missing tool can
  never be reported as a measured miss.

A generative enumerator may well express transformations a rule set cannot. Adding one is
the highest-value next experiment, and the enumerator interface is factored out so that
adding a tool is a single function.

## 7. Files

| File | What it is |
|---|---|
| `run1_generator_recall.py` | The harness. Requires `--structures`; no default input, so it cannot run by accident against a half-filled table. |
| `gold_set_structures_v0.1.csv` | **The real input.** 22 rows, provenance per row. |
| `run1_results_v0.1.csv` | Per-row results at default depth. |
| `run1_results_p3.csv`, `run1_results_p4.csv` | Depth sweep. |
| `gold_set_structures_subset.csv`, `run1_results_subset.csv` | Intermediate 10-row run, superseded. |
| `DEMO_structures_UNVERIFIED.csv`, `DEMO_results_UNVERIFIED.csv` | Original 7 author-entered pipeline-test rows. **Do not cite.** Superseded. |

## 8. Design choices worth knowing

- **Parent-driven rows are excluded from the recall denominator**, not counted as misses.
  Where the toxic species *is* the parent (terfenadine, bosentan, tolcapone, perhexiline),
  recall is not the right question, and the script detects this by comparing InChIKeys
  rather than trusting the `case_type` column.
- **Excluded and unparseable rows are always printed**, never silently dropped —
  otherwise recall quietly inflates as the structure table fills in.
- **Bag size is reported beside recall.** The depth sweep is the reason: recall bought by
  enumerating 657 candidates per compound is not the same result as recall at 82.
