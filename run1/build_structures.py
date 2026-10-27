"""Build run1 structure table. Two provenance routes only:
   A) Wikipedia infobox: SMILES and InChIKey fetched together, so recomputing the
      InChIKey from the SMILES with RDKit verifies transcription.
   B) OPSIN on a chemical name that is RECORDED IN THE GOLD SET (trivial or
      systematic). OPSIN is a deterministic name->structure converter, so the
      only judgement is the name, and the name comes from the curated gold set,
      not from recall.
   Anything outside these two routes is left out and reported as an exclusion.
"""
from py2opsin import py2opsin
from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import MolToInchiKey
import csv, sys
RDLogger.DisableLog("rdApp.*")

WIKI = {
 'acetaminophen':      ('CC(=O)Nc1ccc(O)cc1','RZVAJINKPMORJF-UHFFFAOYSA-N'),
 'cyclophosphamide':   ('O=P1(OCCCN1)N(CCCl)CCCl','CMSMOCZEIVJLDB-UHFFFAOYSA-N'),
 'ifosfamide':         ('O=P1(OCCCN1CCCl)NCCCl','HOMGKSMUEGBAAB-UHFFFAOYSA-N'),
 'terfenadine':        ('OC(c1ccccc1)(c2ccccc2)C4CCN(CCCC(O)c3ccc(cc3)C(C)(C)C)CC4','GUGOEEXESWIERI-UHFFFAOYSA-N'),
 'fexofenadine':       ('O=C(O)C(c1ccc(cc1)C(O)CCCN2CCC(CC2)C(O)(c3ccccc3)c4ccccc4)(C)C','RWTNPBWLLIMQHL-UHFFFAOYSA-N'),
 'capecitabine':       ('FC=1\\C(=N/C(=O)N(C=1)[C@@H]2O[C@@H]([C@@H](O)[C@H]2O)C)\\NC(=O)OCCCCC','GAGWJHPBXLXJQN-UORFTKCHSA-N'),
 'fluorouracil':       ('O=C1NC(=O)NC=C1F','GHASVSINZRGABV-UHFFFAOYSA-N'),
 'irinotecan':         ('O=C7OCC=6C(=O)N2C(\\c1nc5c(c(c1C2)CC)cc(OC(=O)N4CCC(N3CCCCC3)CC4)cc5)=C/C=6[C@@]7(O)CC','UWKQSNNFCGGAFS-XIFFEERXSA-N'),
 'leflunomide':        ('O=C(Nc1ccc(cc1)C(F)(F)F)c2c(onc2)C','VHOGYURTWQBHIL-UHFFFAOYSA-N'),
 'teriflunomide':      ('O=C(Nc1ccc(cc1)C(F)(F)F)C(/C#N)=C(/C)O','UTNUDOFZCWSZMS-YFHOEESVSA-N'),
 'valproic acid':      ('O=C(O)C(CCC)CCC','NIJJYAXOARWZEE-UHFFFAOYSA-N'),
 'isoniazid':          ('C1=CN=CC=C1C(=O)NN','QRXWMOHMRWLFEY-UHFFFAOYSA-N'),
 'carbamazepine':      ('c1ccc2c(c1)C=Cc3ccccc3N2C(=O)N','FFGPTBGBLSHEPO-UHFFFAOYSA-N'),
 'methoxyflurane':     ('ClC(Cl)C(F)(F)OC','RFKMCNOHBTXSMU-UHFFFAOYSA-N'),
 'terbinafine':        ('C(#C\\C=C\\CN(C)Cc2cccc1ccccc12)C(C)(C)C','DOMXUEMWDBAQBQ-WEVVVXLNSA-N'),
 'felbamate':          ('O=C(N)OCC(c1ccccc1)COC(N)=O','WKGXYQFOCVYPAC-UHFFFAOYSA-N'),
 'halothane':          ('BrC(Cl)C(F)(F)F','BCQZXOMGPXTTIC-UHFFFAOYSA-N'),
 'bosentan':           ('CC(C)(C)c1ccc(cc1)S(=O)(=O)Nc2c(c(nc(n2)c3ncccn3)OCCO)Oc4ccccc4OC','GJPICJJJRGTNOD-UHFFFAOYSA-N'),
 'tolcapone':          ('[O-][N+](=O)c2cc(C(=O)c1ccc(cc1)C)cc(O)c2O','MIQPIUSUKVNLNT-UHFFFAOYSA-N'),
 'perhexiline':        ('N3C(CC(C1CCCCC1)C2CCCCC2)CCCC3','CYXKNKQEMFBLER-UHFFFAOYSA-N'),
}
# label -> name as recorded in metabotox_bench_gold_set_v0.2.csv (or its systematic form)
OPSIN = {
 'NAPQI':'N-(4-oxocyclohexa-2,5-dien-1-ylidene)acetamide',
 'acrolein':'prop-2-enal',
 'chloroacetaldehyde':'2-chloroacetaldehyde',
 'methanol':'methanol',
 'formic acid':'methanoic acid',
 'ethylene glycol':'ethane-1,2-diol',
 'oxalic acid':'oxalic acid',
 'glycolic acid':'2-hydroxyacetic acid',
 '4-ene-valproic acid':'2-propylpent-4-enoic acid',
 'hydrazine':'hydrazine',
 'acetylhydrazine':'acetohydrazide',
 'benzene':'benzene',
 'hydroquinone':'benzene-1,4-diol',
 '1,4-benzoquinone':'cyclohexa-2,5-diene-1,4-dione',
 'benzene oxide':'7-oxabicyclo[4.1.0]hepta-2,4-diene',
 'dichloroacetic acid':'2,2-dichloroacetic acid',
 'TBF-A':'(E)-6,6-dimethylhept-2-en-4-ynal',
 'atropaldehyde':'2-phenylprop-2-enal',
 'trifluoroacetyl chloride':'2,2,2-trifluoroacetyl chloride',
}

wiki_ok, prov = {}, {}
print("=== Route A: Wikipedia infobox, InChIKey round-trip check ===")
for k,(smi,key) in sorted(WIKI.items()):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print(f"  {k:24s} SMILES UNPARSEABLE -> excluded"); continue
    calc = MolToInchiKey(m)
    if calc == key:                       st='EXACT'
    elif calc.split('-')[0]==key.split('-')[0]: st='SKELETON-ONLY (stereo differs)'
    else:                                 st='MISMATCH -> EXCLUDED'
    print(f"  {k:24s} {st}")
    if st.startswith('MISMATCH'):
        print(f"      wiki={key}  calc={calc}"); continue
    wiki_ok[k]=Chem.MolToSmiles(m); prov[k]=('wikipedia_infobox', st)

ops_ok={}
print("\n=== Route B: gold-set chemical name -> OPSIN -> SMILES ===")
for label,name in OPSIN.items():
    smi = py2opsin(name)
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None:
        print(f"  {label:24s} OPSIN FAILED on {name!r} -> excluded"); continue
    ops_ok[label]=Chem.MolToSmiles(m); prov[label]=('opsin_from_goldset_name', name)
    print(f"  {label:24s} {ops_ok[label]:40s} {MolToInchiKey(m)}")

def get(k): return wiki_ok.get(k) or ops_ok.get(k)

# (gold-set drug, parent key, toxic species label as in gold set, structure key)
PAIRS = [
 ('Acetaminophen',   'acetaminophen',   'NAPQI (N-acetyl-p-benzoquinone imine)', 'NAPQI'),
 ('Valproic acid',   'valproic acid',   '4-ene-valproic acid',                   '4-ene-valproic acid'),
 ('Ifosfamide',      'ifosfamide',      'chloroacetaldehyde',                    'chloroacetaldehyde'),
 ('Cyclophosphamide','cyclophosphamide','acrolein',                              'acrolein'),
 ('Terbinafine',     'terbinafine',     'TBF-A (E-6,6-dimethylhept-2-en-4-ynal)', 'TBF-A'),
 ('Felbamate',       'felbamate',       'atropaldehyde (2-phenylpropenal)',      'atropaldehyde'),
 ('Halothane',       'halothane',       'trifluoroacetyl chloride',              'trifluoroacetyl chloride'),
 ('Isoniazid',       'isoniazid',       'hydrazine',                             'hydrazine'),
 ('Isoniazid',       'isoniazid',       'acetylhydrazine',                       'acetylhydrazine'),
 ('Capecitabine',    'capecitabine',    '5-fluorouracil (via 5\'-DFUR)',         'fluorouracil'),
 ('Methanol',        'methanol',        'formic acid (via formaldehyde)',        'formic acid'),
 ('Ethylene glycol', 'ethylene glycol', 'glycolic acid',                         'glycolic acid'),
 ('Ethylene glycol', 'ethylene glycol', 'oxalic acid',                           'oxalic acid'),
 ('Methoxyflurane',  'methoxyflurane',  'dichloroacetic acid',                   'dichloroacetic acid'),
 ('Leflunomide',     'leflunomide',     'teriflunomide',                         'teriflunomide'),
 ('Benzene',         'benzene',         'benzene oxide',                         'benzene oxide'),
 ('Benzene',         'benzene',         'hydroquinone',                          'hydroquinone'),
 ('Benzene',         'benzene',         '1,4-benzoquinone',                      '1,4-benzoquinone'),
 # parent-driven: the hazard IS the parent. Recall is not the question; the
 # harness detects these by InChIKey identity and excludes them from the denominator.
 ('Terfenadine',     'terfenadine',     'parent_is_the_hazard',                  'terfenadine'),
 ('Bosentan',        'bosentan',        'parent_is_the_hazard',                  'bosentan'),
 ('Tolcapone',       'tolcapone',       'parent nitrocatechol (uncoupler)',      'tolcapone'),
 ('Perhexiline',     'perhexiline',     'parent accumulates (CYP2D6 PM)',        'perhexiline'),
]
rows, skipped = [], []
for drug,pk,tox_label,tk in PAIRS:
    ps, ts = get(pk), get(tk)
    if not ps or not ts:
        skipped.append((drug,tox_label,'parent' if not ps else 'toxic species')); continue
    pr_p, pr_t = prov.get(pk,('?',''))[0], prov.get(tk,('?',''))[0]
    rows.append(dict(drug_name=drug, parent_smiles=ps, toxic_species_name=tox_label,
                     toxic_species_smiles=ts,
                     notes=f"parent={pr_p}; toxic={pr_t}"))
with open('gold_set_structures_v0.1.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['drug_name','parent_smiles','toxic_species_name','toxic_species_smiles','notes'])
    w.writeheader(); w.writerows(rows)
print(f"\n=== gold_set_structures_v0.1.csv: {len(rows)} rows ===")
for d,t,why in skipped: print(f"  SKIPPED {d} / {t}: no verified {why} structure")
