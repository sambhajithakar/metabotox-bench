# -*- coding: utf-8 -*-
"""Figures for MANUSCRIPT_READY.md. Every number is read from the result files in
../../run1/, ../../run2/ and the gold set, never typed in. Run from this directory.

Figure 1 asserts that all four parent-driven control rows are recovered in every
results file it reads, so a harness regression cannot be plotted silently."""
import csv, collections, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
RUN1 = os.path.join(HERE, "..", "..", "run1")
RUN2 = os.path.join(HERE, "..", "..", "run2")
GOLD = os.path.join(HERE, "..", "..", "metabotox_bench_gold_set_v0.3.csv")

INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8984"
BLUE, ORANGE = "#2a78d6", "#eb6834"     # validated pair: CVD dE 24.7, normal dE 33.6
GRID = "#dcdbd6"

plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 400, "axes.edgecolor": INK3, "axes.linewidth": 0.6,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.labelcolor": INK2, "figure.facecolor": "white", "axes.facecolor": "white",
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

def load(name):
    with open(os.path.join(RUN1, name), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

def stats(rows):
    ev = [r for r in rows if r["status"] == "evaluated"]
    den = [r for r in ev if r.get("parent_is_toxic_species") != "True"]
    hit = sum(1 for r in den if r["sygma_hit_full"] == "True")
    sizes = [int(r["sygma_bag_size"]) for r in den if r.get("sygma_bag_size")]
    return hit, len(den), sum(sizes) / len(sizes)

# ---------------------------------------------------------------- Figure 1
# Coverage against candidate-set COST, two generators, one axis pair.
#
# Cost rather than depth on the x-axis is a deliberate choice, not a stylistic one:
# SyGMa budgets phase-I and phase-II steps separately while an exhaustive rule
# application mixes them, so "depth" is not the same quantity for the two curves.
# Mean candidate-set size is the same quantity for both, and it is the one a user
# pays. This is also why the figure is not two y-axes against depth.
def load2(name):
    with open(os.path.join(RUN2, name), newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))

def stats_gen(rows, hit_col, bag_col):
    ev = [r for r in rows if r["status"] == "evaluated"]
    den = [r for r in ev if r.get("parent_is_toxic_species") != "True"]
    hit = sum(1 for r in den if r[hit_col] == "True")
    sizes = [int(r[bag_col]) for r in den if r.get(bag_col)]
    ctrl = [r for r in ev if r.get("parent_is_toxic_species") == "True"]
    ctrl_ok = sum(1 for r in ctrl if r[hit_col] == "True")
    return hit, len(den), sum(sizes) / len(sizes), ctrl_ok, len(ctrl)

SERIES = [
    ("SyGMa 1.1.0, as distributed", BLUE, "o", [
        ("2 steps", "sygma_p1_1_p2_1.csv", "sygma_hit_full", "sygma_bag_size"),
        ("3 steps\n(default)", "sygma_p1_2_p2_1.csv", "sygma_hit_full", "sygma_bag_size"),
        ("4 steps", "sygma_p1_3_p2_1.csv", "sygma_hit_full", "sygma_bag_size"),
    ]),
    ("GLORYx rule set, exhaustive\n(upper bound on GLORYx)", ORANGE, "s", [
        ("1 step", "run2_results_d1.csv", "gloryx_rules_hit_full", "gloryx_rules_bag_size"),
        ("2 steps", "run2_results_d2.csv", "gloryx_rules_hit_full", "gloryx_rules_bag_size"),
        ("3 steps", "run2_results_d3.csv", "gloryx_rules_hit_full", "gloryx_rules_bag_size"),
    ]),
]

print("Figure 1 points:")
curves = []
for sname, colour, marker, settings in SERIES:
    pts = []
    for label, fn, hc, bc in settings:
        h, n, mean, ck, cn = stats_gen(load2(fn), hc, bc)
        pts.append((label, mean, 100.0 * h / n, h, n))
        print(f"   {sname.splitlines()[0]:28s} {label.replace(chr(10),' '):16s} "
              f"{mean:7.0f} candidates  {h}/{n} = {100.0*h/n:5.1f}%  controls {ck}/{cn}")
        assert ck == cn, f"control row not recovered in {fn} ({ck}/{cn})"
    curves.append((sname, colour, marker, pts))

fig, ax = plt.subplots(figsize=(5.4, 3.3), dpi=400)
# label offsets chosen per point so that no annotation collides with another curve
OFFSETS = {
    "SyGMa 1.1.0, as distributed": [(-2, -13, "top", "center"), (0, -13, "top", "center"), (6, -13, "top", "center")],
    "GLORYx rule set, exhaustive\n(upper bound on GLORYx)":
        [(-10, 1, "bottom", "right"), (0, 12, "bottom", "center"), (-14, 12, "bottom", "center")],
}
for sname, colour, marker, pts in curves:
    xs = [p[1] for p in pts]; ys = [p[2] for p in pts]
    ax.plot(xs, ys, "-", color=colour, lw=2, zorder=2, solid_capstyle="round")
    ax.plot(xs, ys, marker, color=colour, ms=7.5, mec="white", mew=1.8, zorder=3,
            label=sname)
    for (label, m, r, h, n), (dx, dy, va, ha) in zip(pts, OFFSETS[sname]):
        ax.annotate(f"{label}\n{h}/{n} = {r:.1f}%", (m, r), ha=ha, va=va,
                    textcoords="offset points", xytext=(dx, dy), fontsize=7,
                    color=INK2, linespacing=1.3)
ax.set_xscale("log")
ax.set_xticks([10, 20, 50, 100, 200, 500, 1000, 2000])
ax.set_xticklabels(["10", "20", "50", "100", "200", "500", "1000", "2000"])
ax.set_xticks([], minor=True)
ax.set_xlabel("Mean candidate-set size per compound (log scale)")
ax.set_ylabel("Coverage of established\ntoxic species (%)")
ax.set_title("Which rule set you choose is worth ~30 points of coverage at equal cost",
             color=INK, loc="left", pad=8)
ax.set_ylim(0, 100); ax.set_xlim(3.8, 3200)
ax.set_yticks([0, 25, 50, 75, 100])
ax.grid(axis="y", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
ax.legend(loc="lower right", frameon=False, fontsize=7.5, handletextpad=0.5,
          labelspacing=0.8)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "figure1_recall_vs_cost.png"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "figure1_recall_vs_cost.tif"), bbox_inches="tight", dpi=400, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(os.path.join(HERE, "figure1_recall_vs_cost.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
# Per-entry outcome grid: 18 entries x 6 generator settings.
#
# A grid, not a scatter: the claim is about the PATTERN of hits across settings
# (columns that fill in when the rule set changes, rows that stay empty in every
# column), and a pattern across two categorical axes is a matrix. Hit/miss is
# encoded twice, by fill and by shape, so the figure survives greyscale printing
# and colour-vision deficiency.
COLS = [
    ("S2", "sygma_p1_1_p2_1.csv", "sygma_hit_full", "sygma_bag_size"),
    ("S3", "sygma_p1_2_p2_1.csv", "sygma_hit_full", "sygma_bag_size"),
    ("S4", "sygma_p1_3_p2_1.csv", "sygma_hit_full", "sygma_bag_size"),
    ("G1", "run2_results_d1.csv", "gloryx_rules_hit_full", "gloryx_rules_bag_size"),
    ("G2", "run2_results_d2.csv", "gloryx_rules_hit_full", "gloryx_rules_bag_size"),
    ("G3", "run2_results_d3.csv", "gloryx_rules_hit_full", "gloryx_rules_bag_size"),
]
grid, bag_g2 = {}, {}
for cname, fn, hc, bc in COLS:
    for r in load2(fn):
        if r["status"] != "evaluated" or r.get("parent_is_toxic_species") == "True":
            continue
        key = (r["drug_name"], r["toxic_species_name"])
        grid[(key, cname)] = (r[hc] == "True")
        if cname == "G2":
            bag_g2[key] = int(r[bc])
entries = sorted(bag_g2, key=lambda k: (-sum(grid[(k, c[0])] for c in COLS), k[0]))
assert len(entries) == 18, len(entries)
nhit = {c[0]: sum(grid[(k, c[0])] for k in entries) for c in COLS}
print("\nFigure 2 column totals (of 18):", ", ".join(f"{k}={v}" for k, v in nhit.items()))
print("Rows empty in every column:",
      [f"{k[0]} -> {k[1].split(' (')[0]}" for k in entries
       if not any(grid[(k, c[0])] for c in COLS)])

fig, ax = plt.subplots(figsize=(5.9, 4.6), dpi=400)
XS = [0, 1, 2, 3.7, 4.7, 5.7]          # gap separates the two generators
for yy, key in enumerate(entries):
    for xx, (cname, *_rest) in zip(XS, COLS):
        hit = grid[(key, cname)]
        ax.plot([xx], [yy], "s" if hit else "o", ms=8.5, zorder=3,
                color=BLUE if hit else "white",
                mec=BLUE if hit else ORANGE, mew=1.7)
labels = [f"{k[0]} → {k[1].split(' (')[0]}" for k in entries]
ax.set_yticks(range(len(entries)))
ax.set_yticklabels(labels, fontsize=7.5)
# the three rows absent from every setting are the paper's residue: name them in bold
for tick, key in zip(ax.get_yticklabels(), entries):
    if not any(grid[(key, c[0])] for c in COLS):
        tick.set_fontweight("bold"); tick.set_color(INK)
ax.set_ylim(len(entries) - 0.5, -0.5)
ax.set_xticks(XS)
ax.set_xticklabels(["2\nsteps", "3\nsteps", "4\nsteps", "1\nstep", "2\nsteps", "3\nsteps"],
                   fontsize=7.5)
ax.set_xlim(-0.7, 6.4)
ax.tick_params(axis="y", length=0); ax.tick_params(axis="x", length=0)
for sp in ("top", "right", "left", "bottom"): ax.spines[sp].set_visible(False)
ax.grid(axis="y", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
ax.annotate("SyGMa 1.1.0", (1.0, -1.45), xycoords=("data", "data"), ha="center",
            fontsize=8, fontweight="bold", color=INK, annotation_clip=False)
ax.annotate("GLORYx rule set", (4.7, -1.45), xycoords=("data", "data"), ha="center",
            fontsize=8, fontweight="bold", color=INK, annotation_clip=False)
ax.legend(handles=[
    Line2D([], [], marker="s", ls="", ms=8.5, color=BLUE, mec=BLUE, mew=1.7,
           label="species is in the candidate set"),
    Line2D([], [], marker="o", ls="", ms=8.5, color="white", mec=ORANGE, mew=1.7,
           label="absent — unattributable by construction")],
    loc="upper left", frameon=False, fontsize=7.5, handletextpad=0.4,
    bbox_to_anchor=(0.0, -0.085), ncol=2)
fig.suptitle("Nine entries turn on the rule set alone; three turn on nothing",
             x=0.0, ha="left", fontsize=9, color=INK, y=1.02)
fig.savefig(os.path.join(HERE, "figure2_failure_taxonomy.png"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "figure2_failure_taxonomy.tif"), bbox_inches="tight", dpi=400, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(os.path.join(HERE, "figure2_failure_taxonomy.pdf"), bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------- Figure 3
# Cross-organ geometry. Counts are (forming, affected) ORGAN PAIRS, not entries:
# an entry naming two forming or two affected organs contributes several pairs.
rows = list(csv.DictReader(open(GOLD, encoding="utf-8")))
same, cross, pairs = collections.Counter(), collections.Counter(), 0
for r in rows:
    for f in [x.strip() for x in r["forming_organ"].split(";") if x.strip()]:
        for a in [x.strip() for x in r["affected_organ"].split(";") if x.strip()]:
            if a == "none": continue
            pairs += 1
            (same if f == a else cross)[a] += 1
affected = sorted(set(same) | set(cross), key=lambda a: -(same[a] + cross[a]))
print(f"\nFigure 3: {pairs} (forming, affected) organ pairs across {len(rows)} entries; "
      f"{sum(cross.values())} cross-organ, {sum(same.values())} same-organ")

fig, ax = plt.subplots(figsize=(4.6, 4.3), dpi=400)
ys = range(len(affected))
s = [same[a] for a in affected]; c = [cross[a] for a in affected]
ax.barh(ys, s, color=BLUE, height=0.62, label="Formed in the organ it injures", zorder=3)
ax.barh(ys, c, left=s, color=ORANGE, height=0.62, hatch="///", edgecolor="white",
        linewidth=0.9, label="Formed in a different organ", zorder=3)
for i, a in enumerate(affected):
    t = same[a] + cross[a]
    ax.annotate(f"{t}", (t, i), textcoords="offset points", xytext=(5, 0),
                fontsize=7, color=INK3, va="center")
ax.set_yticks(list(ys)); ax.set_yticklabels([a.replace("_", " ") for a in affected], fontsize=7.5)
ax.invert_yaxis()
ax.set_xlabel("(forming organ, affected organ) pairs in MetaboTox-Bench v0.3")
ax.set_title("Outside the liver, injury is almost always cross-organ", color=INK, loc="left", pad=8)
ax.grid(axis="x", color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.legend(loc="lower right", frameon=False, fontsize=7.5)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figure3_cross_organ.png"), bbox_inches="tight")
fig.savefig(os.path.join(HERE, "figure3_cross_organ.tif"), bbox_inches="tight", dpi=400, pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(os.path.join(HERE, "figure3_cross_organ.pdf"), bbox_inches="tight"); plt.close(fig)

print("\nWrote figure1/2/3 as .png (400 dpi), .tif (400 dpi, LZW) and .pdf (vector) in", HERE)
