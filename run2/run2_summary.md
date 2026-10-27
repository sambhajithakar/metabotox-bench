# Run 2 — GLORYx rule-set coverage vs SyGMa, on one structure table

Executed 13 October 2026. Harness: `run2_gloryx_rules_recall.py` (GLORYx rules) and
`../run1/run1_generator_recall.py` (SyGMa). Both use the same InChIKey matching, the
same exclusion handling and the same recall denominator. Structure table:
`../run1/gold_set_structures_v0.1.csv` (22 rows: 18 metabolite-driven, 4 parent-driven
controls).

## What was and was not run

**GLORYx itself was not run, and cannot be from public artifacts.** GLORYx = FAME 3
site-of-metabolism prediction + reaction rules. The repository
(`github.com/christinadebruynkops/gloryx`, commit `a1f4a824`) states that the FAME 3
models are not bundled and must be obtained from the FAME 3 authors. The reaction rules
*are* published, in `reaction_rules/gloryx_reactionrules.csv` (copied here verbatim,
md5 `be856a6a4a58db8a5ebfc8df43158ad1`).

What was measured is the **expressive ceiling of that rule set**: all 260 rules applied
exhaustively, with no site-of-metabolism filtering, no priority filtering and no scoring.
Because FAME 3 can only select among products the rules generate, this is a generous
**upper bound** on GLORYx's own coverage. A miss here is a miss GLORYx cannot avoid; a
hit here is not a prediction GLORYx would necessarily make.

All 260 SMIRKS compiled under RDKit 2026.09.1 (0 parse failures). Rule provenance as
distributed: 178 from SyGMa, 73 from GLORY, 9 from the GLORYx authors; 206 `common`,
54 `uncommon`; 237 phase 1, 23 phase 2 — all applied.

## Headline result

Coverage of the established toxic species, n = 18, by generator and step budget:

| Generator | Rule applications | Coverage (full InChIKey) | Coverage (skeleton) | Mean candidate-set size |
|---|---|---|---|---|
| SyGMa (phase1=1, phase2=1) | 2 | 4/18 = 22.2% | 4/18 = 22.2% | 18 |
| SyGMa (phase1=2, phase2=1) — documented default | 3 | 7/18 = 38.9% | 7/18 = 38.9% | 82 |
| SyGMa (phase1=3, phase2=1) | 4 | 8/18 = 44.4% | 8/18 = 44.4% | 263 |
| GLORYx rule set, exhaustive | 1 | 7/18 = 38.9% | 7/18 = 38.9% | 14 |
| GLORYx rule set, exhaustive | 2 | 13/18 = 72.2% | 13/18 = 72.2% | 170 |
| GLORYx rule set, exhaustive | 3 | 15/18 = 83.3% | 15/18 = 83.3% | 1992 |

Full-key and skeleton coverage are identical at every setting and for both generators, so
every miss is a structural absence, not a stereochemical near miss. No row reached the
frontier cap at depth 1–3, so none of those figures is a cap artefact. All four
parent-driven controls are recovered by the GLORYx rule set at depth 2, as by SyGMa,
confirming the matching procedure.

Two readings, both of which belong in the paper:

* **At a matched step budget**, the GLORYx rule set covers 72.2% at two applications
  against SyGMa's 22.2%.
* **At a matched candidate-set cost** — the quantity that actually governs downstream
  screening burden — the GLORYx rule set covers 72.2% with 170 candidates per compound,
  against SyGMa's 44.4% with 263.

## The three residual misses are depth-invariant

| Entry | d3 | d4 | d5 | candidate set at d5 |
|---|---|---|---|---|
| Methoxyflurane → dichloroacetic acid | miss | miss | miss | 29,972 (no cap) |
| Leflunomide → teriflunomide | miss | miss | miss | 60,003 (**cap reached**, so a lower bound) |
| Benzene → benzene oxide | miss | miss | miss | 6,937 (no cap) |

Methoxyflurane and benzene are therefore absences established with certainty: thirty
thousand and seven thousand distinct products respectively, and the species is in neither.
Leflunomide was cap-free and still absent at depth 4 (17,373 products).

What each failure is:

* **Methoxyflurane → dichloroacetic acid.** Oxidative dehalogenation of a polyhalogenated
  ether. Compare halothane → trifluoroacetyl chloride, which the GLORYx rules *do* reach
  at depth 1 and SyGMa never does, so dehalogenation is partly encoded and partly not.
* **Leflunomide → teriflunomide.** Isoxazole ring-opening. Rule sets express ring
  oxidation and ring hydroxylation but not the reductive ring cleavage that produces the
  active, and hepatotoxic, open-chain enol-nitrile.
* **Benzene → benzene oxide.** The most instructive of the three. Phenol, hydroquinone and
  1,4-benzoquinone are all reached by depth 2, so benzene's downstream metabolites are
  well covered — but the arene oxide itself, which is the electrophile, is not. The rules
  encode the **net transformation** (C–H to C–OH) and step over the intermediate.

That last pattern is the general lesson and it is not specific to benzene: a rule set
written to reproduce observed, isolable metabolites has no reason to emit a transient
intermediate, and the transient intermediate is frequently the toxicologically
responsible species.

## Files

| File | What it is |
|---|---|
| `gloryx_reactionrules.csv` | the 260 published GLORYx rules, copied verbatim |
| `GLORYX_SOURCE_COMMIT.txt` | source repository and commit the rules came from |
| `run2_gloryx_rules_recall.py` | the coverage harness |
| `run2_results_d1.csv`, `_d2.csv`, `_d3.csv` | per-row results at depth 1, 2, 3 |
| `residual_misses.csv` | the three-row subset used for the depth-4/5 test |
| `run2_residual_d4.csv`, `_d5.csv` | per-row results for that subset |
| `run2_d3.log`, `run2_residual_d4_d5.log` | console logs including self-check output |
| `sygma_p1_1_p2_1.csv`, `sygma_p1_2_p2_1.csv`, `sygma_p1_3_p2_1.csv` | the matched SyGMa re-runs |

## Two harness decisions worth recording

**The parent is retained in the candidate set, in both harnesses.** SyGMa's `tree.to_list()`
includes the root molecule, so the four parent-driven control rows — where the hazard *is*
the administered compound — work as a positive control on the InChIKey matching. An earlier
version of the GLORYx harness subtracted the parent, which made those same four rows fail by
construction and destroyed the control. Re-running with the parent retained changed no
outcome among the 18 metabolite-driven rows (the toxic species is non-identical to the parent
there by the exclusion rule, so the parent can never create a false hit) and raised each
candidate-set size by exactly one. The figures above are the parent-retained run; the
pre-fix run gave 38.9% / 72.2% / 83.3% as well, with means 13 / 169 / 1991.

## One bug worth recording

The first version of the harness reported 5.6% coverage at depth 1 with a mean candidate
set of 4. That was wrong. GLORYx's SMIRKS are written for CDK and name hydrogen as an
explicit atom in the reactant template (`[cH1:5]([H])`), so on an implicit-hydrogen
molecule almost every rule silently fails to match. Adding `Chem.AddHs` before
`RunReactants` took paracetamol from 1 product to 13, NAPQI among them, and depth-1
coverage from 5.6% to 38.9%.

The harness now refuses to print any coverage number unless a self-check passes:
paracetamol must yield at least 10 products in one step and NAPQI must be among them.
A silent collapse in mean candidate-set size is the symptom to watch for if this
regression ever returns.
