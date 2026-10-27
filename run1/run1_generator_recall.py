#!/usr/bin/env python
"""
Run 1 — metabolite generator recall ceiling for MetaboTox-Bench.

Question: of the established toxic species in the gold set, how many does a
metabolite enumerator actually propose? If a generator never proposes NAPQI,
no downstream model can attribute acetaminophen hepatotoxicity to it. This is
a hard upper bound on every attribution metric in the paper, which is why it
is Run 1 and why it is reported before any model result.

Status: SyGMa is wired up and verified working. BioTransformer (Java) and
MetaPredictor are stubs — see run1_README.md.

Input CSV (--structures) must have columns:
    drug_name,parent_smiles,toxic_species_name,toxic_species_smiles,notes
One row per (drug, toxic species). A drug with two toxic species gets two rows.
Rows with an empty SMILES field are counted as `missing_structure` and excluded
from recall, with the exclusion reported — never silently dropped.

Matching is by InChIKey. Two modes are reported because they answer different
questions:
    full      - full InChIKey. Exact structure including stereochemistry.
    skeleton  - first 14-character block only. Ignores stereochemistry and
                protonation. Enumerators routinely get connectivity right and
                stereo wrong, so skeleton recall is the fairer ceiling and
                full recall is the strict one. Report both.

Usage:
    python run1_generator_recall.py --structures gold_set_structures.csv \
        --phases 2 --out run1_results.csv

Deliberately has no default input path: it must not be runnable by accident
against a half-filled structure table.
"""

import argparse
import csv
import sys
from collections import Counter

try:
    from rdkit import Chem
    from rdkit import RDLogger
    from rdkit.Chem.inchi import MolToInchiKey, INCHI_AVAILABLE
except ImportError:
    sys.exit("RDKit not importable. See run1_README.md for the install that works here.")

if not INCHI_AVAILABLE:
    sys.exit("RDKit lacks InChI support; matching is by InChIKey so this cannot proceed.")

RDLogger.DisableLog("rdApp.*")


def inchikey(smiles):
    """InChIKey for a SMILES string, or None if it will not parse."""
    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None
    try:
        return MolToInchiKey(mol)
    except Exception:
        return None


def skeleton(key):
    return key.split("-")[0] if key else None


# ---------------------------------------------------------------- enumerators

def enumerate_sygma(parent_smiles, phase1=2, phase2=1):
    """SyGMa phase I/II enumeration. Returns {InChIKey: score}.

    Default 2 phase-I steps and 1 phase-II step follows SyGMa's own documented
    example. Raising these explodes the bag size, so if you change them, report
    mean bag size alongside recall — recall bought with a 10,000-member bag is
    not a usable result.
    """
    import sygma

    scenario = sygma.Scenario([
        [sygma.ruleset["phase1"], phase1],
        [sygma.ruleset["phase2"], phase2],
    ])
    mol = Chem.MolFromSmiles(parent_smiles)
    if mol is None:
        return None  # signals unparseable parent, distinct from empty result
    tree = scenario.run(mol)
    # calc_scores() must precede to_list(): to_list sorts on the score field,
    # which is None until scores are computed, and raises TypeError otherwise.
    tree.calc_scores()
    out = {}
    for row in tree.to_list():
        m = row.get("SyGMa_metabolite")
        if m is None:
            continue
        try:
            k = MolToInchiKey(m)
        except Exception:
            continue
        score = float(row.get("SyGMa_score") or 0.0)
        if k and score > out.get(k, -1.0):
            out[k] = score
    return out


def enumerate_biotransformer(parent_smiles, **kw):
    """NOT IMPLEMENTED. Needs the BioTransformer 3 jar; java is present here.

    Do not substitute a different tool and report it under this name.
    """
    raise NotImplementedError("BioTransformer jar not present; see run1_README.md")


def enumerate_metapredictor(parent_smiles, **kw):
    """NOT IMPLEMENTED. Needs model weights not available in this container."""
    raise NotImplementedError("MetaPredictor weights not present; see run1_README.md")


ENUMERATORS = {
    "sygma": enumerate_sygma,
    "biotransformer": enumerate_biotransformer,
    "metapredictor": enumerate_metapredictor,
}


# ------------------------------------------------------------------- main run

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structures", required=True,
                    help="CSV: drug_name,parent_smiles,toxic_species_name,toxic_species_smiles,notes")
    ap.add_argument("--out", required=True, help="per-row results CSV")
    ap.add_argument("--enumerators", default="sygma",
                    help="comma-separated subset of: " + ",".join(ENUMERATORS))
    ap.add_argument("--phase1", type=int, default=2)
    ap.add_argument("--phase2", type=int, default=1)
    args = ap.parse_args()

    names = [e.strip() for e in args.enumerators.split(",") if e.strip()]
    for n in names:
        if n not in ENUMERATORS:
            sys.exit(f"unknown enumerator {n!r}; choose from {list(ENUMERATORS)}")

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
        # A toxic species identical to the parent is a parent-driven case
        # (bosentan, tolcapone). Recall is not the right question there.
        rec["parent_is_toxic_species"] = str(parent_key == tox_key)

        for n in names:
            try:
                bag = ENUMERATORS[n](parent, phase1=args.phase1, phase2=args.phase2)
            except NotImplementedError as e:
                rec[f"{n}_hit_full"] = "not_implemented"
                rec[f"{n}_hit_skeleton"] = "not_implemented"
                rec[f"{n}_bag_size"] = ""
                rec[f"{n}_note"] = str(e)
                continue
            except Exception as e:  # a crash on one molecule must not lose the run
                rec[f"{n}_hit_full"] = "error"
                rec[f"{n}_hit_skeleton"] = "error"
                rec[f"{n}_bag_size"] = ""
                rec[f"{n}_note"] = f"{type(e).__name__}: {e}"[:200]
                continue
            if bag is None:
                rec[f"{n}_hit_full"] = "error"
                rec[f"{n}_hit_skeleton"] = "error"
                rec[f"{n}_bag_size"] = ""
                rec[f"{n}_note"] = "parent unparseable by enumerator"
                continue
            skels = {skeleton(k) for k in bag}
            rec[f"{n}_hit_full"] = str(tox_key in bag)
            rec[f"{n}_hit_skeleton"] = str(skeleton(tox_key) in skels)
            rec[f"{n}_bag_size"] = str(len(bag))
            rec[f"{n}_note"] = ""
        results.append(rec)

    cols = ["drug_name", "toxic_species_name", "status", "toxic_inchikey",
            "parent_is_toxic_species"]
    for n in names:
        cols += [f"{n}_hit_full", f"{n}_hit_skeleton", f"{n}_bag_size", f"{n}_note"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(results)

    # -------- summary. Excluded rows are always printed, never hidden.
    ev = [r for r in results if r["status"] == "evaluated"]
    excl = Counter(r["status"] for r in results if r["status"] != "evaluated")
    print(f"\nRows in table:        {len(results)}")
    print(f"Evaluated:            {len(ev)}")
    for k, v in excl.items():
        print(f"EXCLUDED [{k}]: {v}")
    pdriven = [r for r in ev if r.get("parent_is_toxic_species") == "True"]
    if pdriven:
        print(f"Parent-driven rows (recall not applicable): {len(pdriven)} "
              f"-> {', '.join(r['drug_name'] for r in pdriven)}")

    denom = [r for r in ev if r.get("parent_is_toxic_species") != "True"]
    print(f"Recall denominator:   {len(denom)}")
    for n in names:
        full = sum(1 for r in denom if r.get(f"{n}_hit_full") == "True")
        skel = sum(1 for r in denom if r.get(f"{n}_hit_skeleton") == "True")
        errs = sum(1 for r in denom if r.get(f"{n}_hit_full") in ("error", "not_implemented"))
        sizes = [int(r[f"{n}_bag_size"]) for r in denom if r.get(f"{n}_bag_size")]
        mean = sum(sizes) / len(sizes) if sizes else 0
        if denom:
            print(f"\n  {n}:")
            print(f"    recall (full InChIKey):     {full}/{len(denom)} = {full/len(denom):.1%}")
            print(f"    recall (skeleton, no stereo): {skel}/{len(denom)} = {skel/len(denom):.1%}")
            print(f"    mean bag size:              {mean:.0f}")
            if errs:
                print(f"    errors / not implemented:   {errs}")
        misses = [r["drug_name"] for r in denom if r.get(f"{n}_hit_skeleton") == "False"]
        if misses:
            print(f"    MISSED (unattributable by construction): {', '.join(misses)}")
    print("\nMissed species are the ceiling. Report them in Table 1 and state in every")
    print("attribution table whether they are included or excluded.\n")


if __name__ == "__main__":
    main()
