#!/usr/bin/env python
"""
Run 2 — coverage ceiling of the GLORYx reaction rule set.

WHAT THIS IS, AND WHAT IT IS NOT
--------------------------------
This script does NOT run GLORYx. GLORYx = FAME 3 site-of-metabolism prediction
+ reaction rules, and the FAME 3 models are explicitly not distributed with the
GLORYx source ("the models need to be obtained from the authors of FAME 3"),
so GLORYx cannot be executed from public artifacts alone. That is reported as a
reproducibility finding, not worked around.

What this script DOES measure is the *expressive ceiling* of GLORYx's published
reaction rule set: it applies all of GLORYx's reaction rules exhaustively, with
no site-of-metabolism filtering and no scoring, and asks whether the resulting
product space contains each established toxic species.

Why that is the right quantity for the paper's claim. In GLORYx, FAME 3 selects
and ranks among the products that the rule set can generate; it cannot invent a
product no rule expresses. Therefore

    recall(GLORYx) <= coverage(GLORYx rule set, exhaustive)

and the number this script produces is a GENEROUS UPPER BOUND on what GLORYx
itself could achieve on this table. A miss here is a miss GLORYx cannot avoid.
A hit here is not a prediction GLORYx would necessarily make.

Every deliberate choice below is biased in the rule set's favour, so that the
bound stays an upper bound:
  * no SoM filtering, no priority filtering: `uncommon` rules are applied too
  * both phase 1 and phase 2 rule subsets applied at every step
  * all product fragments retained, including small ones
  * matching on the 14-character InChIKey skeleton as well as the full key

Matching, exclusion handling and the recall denominator are identical to
run1_generator_recall.py, so SyGMa and the GLORYx rule set are compared under
one routine on one structure table.

Usage:
    python run2_gloryx_rules_recall.py \
        --structures ../run1/gold_set_structures_v0.1.csv \
        --rules gloryx_reactionrules.csv \
        --depth 2 --out run2_results_d2.csv
"""

import argparse
import csv
import sys
import time
from collections import Counter

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    from rdkit.Chem.inchi import MolToInchiKey, INCHI_AVAILABLE
except ImportError:
    sys.exit("RDKit not importable. See ../run1/run1_README.md for the install that works here.")

if not INCHI_AVAILABLE:
    sys.exit("RDKit lacks InChI support; matching is by InChIKey so this cannot proceed.")

RDLogger.DisableLog("rdApp.*")

# A frontier cap is needed because exhaustive application of 260 rules explodes.
# If the cap is reached for a molecule the result for that molecule is a LOWER
# bound on rule-set coverage, which is recorded per row so it is never silently
# conflated with a genuine absence.
DEFAULT_FRONTIER_CAP = 20000


def inchikey_from_mol(mol):
    try:
        return MolToInchiKey(mol)
    except Exception:
        return None


def inchikey(smiles):
    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    return inchikey_from_mol(mol)


def skeleton(key):
    return key.split("-")[0] if key else None


def load_rules(path):
    """Returns (reactions, stats). A SMIRKS RDKit cannot compile is counted and
    named, never skipped quietly: an uncounted parse failure would understate
    the rule set and so overstate this paper's claim."""
    reactions, failed = [], []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            smirks = (row.get("SMIRKS") or "").strip()
            name = (row.get("Reaction name") or "").strip()
            if not smirks:
                failed.append((name, "empty SMIRKS"))
                continue
            try:
                rxn = AllChem.ReactionFromSmarts(smirks)
                if rxn is None:
                    failed.append((name, "ReactionFromSmarts returned None"))
                    continue
                rxn.Initialize()
            except Exception as e:
                failed.append((name, f"{type(e).__name__}: {e}"[:120]))
                continue
            if rxn.GetNumReactantTemplates() != 1:
                # Every GLORYx rule is unimolecular on the xenobiotic; a
                # multi-reactant template would need a cofactor we do not supply.
                failed.append((name, f"{rxn.GetNumReactantTemplates()} reactant templates"))
                continue
            reactions.append((name, row.get("Name of rule subset", ""),
                              row.get("Priority level", ""), rxn))
    return reactions, failed


def apply_rules_once(mol, reactions):
    """All single-step products of `mol` under every rule. Fragments split and
    kept separately, which is the generous reading.

    EXPLICIT HYDROGENS ARE LOAD-BEARING. GLORYx's SMIRKS are written for CDK
    and name hydrogen as a real atom in the reactant template, e.g.
        [#6:1]~[a:2]1[a:3][a:4][cH1:5]([H])[a:6][a:7]1>>...(O)...
    On an implicit-H molecule the `([H])` branch cannot match, so almost every
    rule silently fails and the measured coverage collapses. Measured here: on
    paracetamol, 1 product without AddHs versus 13 with it, and NAPQI appears
    only in the second case. Any future variant of this harness must keep the
    AddHs, and a sudden drop in mean candidate-set size is the symptom to look
    for if it is lost.

    Products are round-tripped through canonical SMILES so that the next
    generation starts from consistently valenced molecules rather than from
    RDKit's raw reaction output.
    """
    out = {}
    mol_h = Chem.AddHs(mol)
    for _name, _subset, _prio, rxn in reactions:
        try:
            product_sets = rxn.RunReactants((mol_h,))
        except Exception:
            continue
        for products in product_sets:
            for p in products:
                try:
                    p.UpdatePropertyCache(strict=False)
                    Chem.SanitizeMol(p)
                    clean = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.RemoveHs(p)))
                except Exception:
                    continue
                if clean is None:
                    continue
                try:
                    frags = Chem.GetMolFrags(clean, asMols=True)
                except Exception:
                    continue
                for frag in frags:
                    k = inchikey_from_mol(frag)
                    if k and k not in out:
                        out[k] = frag
    return out


def enumerate_gloryx_rules(parent_smiles, reactions, depth, frontier_cap):
    """Exhaustive product space to `depth` rule applications.

    Returns (bag_keys, capped, n_generations_run).

    The parent is RETAINED in the returned bag. That is deliberate and it is what
    makes the two harnesses comparable: SyGMa's `tree.to_list()` includes the root
    molecule, so the four parent-driven control rows — where the hazard *is* the
    administered compound — are a positive control on the InChIKey matching for
    SyGMa. Excluding the parent here would have made those same four rows fail by
    construction and destroyed the control, while leaving the 18 metabolite-driven
    rows untouched (on those the toxic species is non-identical to the parent by
    the exclusion rule in §4.1, so the parent can never create a false hit).
    Verified empirically: including the parent changes no outcome among the 18 and
    raises each candidate-set size by exactly one.
    """
    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        return None, False, 0
    parent_key = inchikey_from_mol(mol)
    seen = {parent_key: mol} if parent_key else {}
    frontier = {parent_key: mol} if parent_key else {}
    capped = False
    gens = 0
    for _ in range(depth):
        new = {}
        for _k, m in frontier.items():
            for k2, m2 in apply_rules_once(m, reactions).items():
                if k2 not in seen and k2 not in new:
                    new[k2] = m2
            if len(seen) + len(new) > frontier_cap:
                capped = True
                break
        if not new:
            break
        seen.update(new)
        frontier = new
        gens += 1
        if capped:
            break
    return set(seen), capped, gens


#: Paracetamol, and the InChIKey skeleton of NAPQI. Used as a tripwire: the
#: rule set must turn the first into the second in one step. This caught a real
#: bug (missing AddHs, see apply_rules_once) that had silently cut measured
#: coverage from 44% to 6%, which would have been published as a result.
SELFCHECK_PARENT = "CC(=O)Nc1ccc(O)cc1"
SELFCHECK_PRODUCT_SKELETON = skeleton(inchikey("CC(=O)N=C1C=CC(=O)C=C1"))
SELFCHECK_MIN_PRODUCTS = 10


def selfcheck(reactions):
    """Abort rather than report a number produced by a broken rule application."""
    mol = Chem.MolFromSmiles(SELFCHECK_PARENT)
    bag = apply_rules_once(mol, reactions)
    skels = {skeleton(k) for k in bag}
    ok_size = len(bag) >= SELFCHECK_MIN_PRODUCTS
    ok_napqi = SELFCHECK_PRODUCT_SKELETON in skels
    print(f"self-check: paracetamol -> {len(bag)} products in one step "
          f"(need >={SELFCHECK_MIN_PRODUCTS}), NAPQI present: {ok_napqi}")
    if not (ok_size and ok_napqi):
        sys.exit("SELF-CHECK FAILED: rule application is broken (explicit "
                 "hydrogens? rule file?). Refusing to report coverage numbers.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", required=True)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--depth", type=int, default=2,
                    help="number of sequential rule applications")
    ap.add_argument("--frontier-cap", type=int, default=DEFAULT_FRONTIER_CAP)
    args = ap.parse_args()

    reactions, failed = load_rules(args.rules)
    print(f"GLORYx rules compiled: {len(reactions)}")
    if failed:
        print(f"GLORYx rules NOT compiled: {len(failed)}")
        for name, why in failed:
            print(f"    - {name}: {why}")
    else:
        print("GLORYx rules NOT compiled: 0")
    selfcheck(reactions)

    with open(args.structures, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit("structure table is empty")

    results = []
    for r in rows:
        drug = (r.get("drug_name") or "").strip()
        parent = (r.get("parent_smiles") or "").strip()
        tox_name = (r.get("toxic_species_name") or "").strip()
        tox_smi = (r.get("toxic_species_smiles") or "").strip()

        rec = {"drug_name": drug, "toxic_species_name": tox_name}
        tox_key = inchikey(tox_smi)
        parent_key = inchikey(parent)

        if not parent or parent_key is None:
            rec["status"] = "missing_or_bad_parent_structure"
            results.append(rec)
            continue
        if not tox_smi or tox_key is None:
            rec["status"] = "missing_or_bad_toxic_structure"
            results.append(rec)
            continue

        rec["status"] = "evaluated"
        rec["toxic_inchikey"] = tox_key
        rec["parent_is_toxic_species"] = str(parent_key == tox_key)

        t0 = time.time()
        bag, capped, gens = enumerate_gloryx_rules(
            parent, reactions, args.depth, args.frontier_cap)
        dt = time.time() - t0
        if bag is None:
            rec["gloryx_rules_hit_full"] = "error"
            rec["gloryx_rules_hit_skeleton"] = "error"
            rec["gloryx_rules_bag_size"] = ""
            rec["gloryx_rules_note"] = "parent unparseable"
        else:
            skels = {skeleton(k) for k in bag}
            rec["gloryx_rules_hit_full"] = str(tox_key in bag)
            rec["gloryx_rules_hit_skeleton"] = str(skeleton(tox_key) in skels)
            rec["gloryx_rules_bag_size"] = str(len(bag))
            rec["gloryx_rules_note"] = (
                f"frontier_cap_reached generations={gens}" if capped else "")
        rec["seconds"] = f"{dt:.1f}"
        print(f"  {drug:22s} {tox_name[:34]:34s} "
              f"bag={rec.get('gloryx_rules_bag_size',''):>6s} "
              f"hit={rec.get('gloryx_rules_hit_skeleton','')} "
              f"{rec['gloryx_rules_note']} [{dt:.0f}s]", flush=True)
        results.append(rec)

    cols = ["drug_name", "toxic_species_name", "status", "toxic_inchikey",
            "parent_is_toxic_species", "gloryx_rules_hit_full",
            "gloryx_rules_hit_skeleton", "gloryx_rules_bag_size",
            "gloryx_rules_note", "seconds"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    ev = [r for r in results if r["status"] == "evaluated"]
    excl = Counter(r["status"] for r in results if r["status"] != "evaluated")
    print(f"\nDepth:                {args.depth}")
    print(f"Rows in table:        {len(results)}")
    print(f"Evaluated:            {len(ev)}")
    for k, v in excl.items():
        print(f"EXCLUDED [{k}]: {v}")
    pdriven = [r for r in ev if r.get("parent_is_toxic_species") == "True"]
    if pdriven:
        print(f"Parent-driven rows (recall not applicable): {len(pdriven)} "
              f"-> {', '.join(r['drug_name'] for r in pdriven)}")

    denom = [r for r in ev if r.get("parent_is_toxic_species") != "True"]
    full = sum(1 for r in denom if r.get("gloryx_rules_hit_full") == "True")
    skel = sum(1 for r in denom if r.get("gloryx_rules_hit_skeleton") == "True")
    sizes = [int(r["gloryx_rules_bag_size"]) for r in denom if r.get("gloryx_rules_bag_size")]
    mean = sum(sizes) / len(sizes) if sizes else 0
    ncap = sum(1 for r in denom if "frontier_cap_reached" in (r.get("gloryx_rules_note") or ""))
    print(f"Recall denominator:   {len(denom)}")
    print(f"\n  gloryx_rules (exhaustive, UPPER BOUND on GLORYx):")
    print(f"    coverage (full InChIKey):     {full}/{len(denom)} = {full/len(denom):.1%}")
    print(f"    coverage (skeleton):          {skel}/{len(denom)} = {skel/len(denom):.1%}")
    print(f"    mean candidate-set size:      {mean:.0f}")
    print(f"    rows hitting frontier cap:    {ncap}"
          + ("  (those rows are lower bounds)" if ncap else ""))
    misses = [r["drug_name"] + "/" + r["toxic_species_name"][:24]
              for r in denom if r.get("gloryx_rules_hit_skeleton") == "False"]
    if misses:
        print("    NOT EXPRESSIBLE by any GLORYx rule within this depth:")
        for m in misses:
            print(f"      - {m}")


if __name__ == "__main__":
    main()
