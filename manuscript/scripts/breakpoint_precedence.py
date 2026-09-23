"""Which caller sets the breakpoints of a consensus call.

Where used
----------
Supplementary Information (`manuscript/supplementary.typ`):
    Supplementary Note S3                   every number in the note
    Supplementary Table Breakpoint Precedence  precedence and span, all coverages
    Supplementary Figure Breakpoint Precedence A, pairwise outward extension;
                                               B, sole boundary precedence by class
Discussion -> the union-merging paragraph, which states the note's conclusion
without numbers and cites it.

A union merge gives a component the minimum start and maximum end of its member
calls, so in a component carrying two or more callers the coordinates come from
whichever caller reached furthest in each direction. Two quantities settle which
caller that is, and they are not the same question:

  * *pairwise outward extension*: for two callers in one component, how far
    beyond the other each one reaches at each edge. This is the permissiveness
    of one caller relative to another, in base pairs and in sign.
  * *boundary precedence*: over all components the caller appears in, how often
    it is the one holding the union's start or end. A caller can take precedence
    without being much wider, because CNVpytor and GATK-gCNV are both quantized
    to the 1 kb bin and therefore tie with each other exactly, while Delly's
    split-read breakpoints almost never land on the grid and so break the tie in
    one direction or the other.

Components are the 30x consensus merge at the adopted 0.5 reciprocal overlap,
restricted to those carrying at least two callers; the other three coverages are
reported in the table for robustness. A caller contributing several member calls
to one component is summarized by its own minimum start and maximum end, which
is what that caller alone would have contributed to the union.

    pixi run python manuscript/scripts/breakpoint_precedence.py
"""

import glob
import itertools
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

from consensuscnv.callsets import collect_callsets, merge_components
from consensuscnv.callsets.registry import SOURCES, SVTYPES, seed_chromosomes
from consensuscnv.utils import read_genome_file

ROOT = Path("/lab01/Projects/Lionel_Projects/blendedCNV_pipeline")

# Chromosome ids have to order the genome rather than the order the BEDs
# happen to be read in, so the registry is seeded before any CallSet is built.
seed_chromosomes(read_genome_file(ROOT / "src" / "consensuscnv" / "templates" / "genome_primary_hg38.txt"))
DEST = ROOT / "results" / "breakpoint_precedence"
TABLES = ROOT / "results" / "manuscript"

COVERAGES = ("30x", "6x", "4x", "2x")
CALLERS = ("cnvpytor", "delly", "gatk")
LABELS = {"cnvpytor": "CNVpytor", "delly": "Delly", "gatk": "GATK-gCNV"}
CONSENSUS_THRESHOLD = 0.5
FOCUS = "30x"

# Okabe-Ito for the callers, as in Figure 4 and Supplementary Figure Size Floor.
COLORS = {"CNVpytor": "#0072B2", "Delly": "#D55E00", "GATK-gCNV": "#009E73"}

mpl.rcParams.update({
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "font.family": "DejaVu Sans", "figure.dpi": 150,
})


def bed_paths(coverage: str) -> list[str]:
    """Every per-sample BED for the three callers at one coverage.

    Pinned to the caller directories rather than a wildcard: consensus output
    lands beside them, and a wildcard would read merged calls back in as an
    extra source.
    """
    root = ROOT / "out" / f"{coverage}_Coverage"
    return [path for caller in CALLERS for path in sorted(glob.glob(str(root / caller / "*.bed")))]


def per_caller_extents(coverage: str):
    """Each component's per-caller extents, for components carrying >= 2 callers.

    Returns `(lo, hi, present, svtype)`, the first three of shape
    (components, callers): the minimum start and maximum end each caller
    contributed, and whether it contributed at all.
    """
    calls = collect_callsets(bed_paths(coverage))
    merged = merge_components(calls, min_reciprocal_overlap=CONSENSUS_THRESHOLD)

    # A leaf call carries exactly one source bit, so its position is the index.
    source_idx = np.log2(calls.source_bits).astype(np.int64)
    shape = (len(merged), len(SOURCES.names))
    lo = np.full(shape, np.iinfo(np.int64).max)
    hi = np.full(shape, np.iinfo(np.int64).min)
    for index in range(shape[1]):
        rows = source_idx == index
        np.minimum.at(lo[:, index], merged.labels[rows], calls.starts[rows])
        np.maximum.at(hi[:, index], merged.labels[rows], calls.ends[rows])

    present = lo != np.iinfo(np.int64).max
    multi = present.sum(axis=1) >= 2
    svtype = np.array(SVTYPES.names)[merged.svtype_idx]
    return lo[multi], hi[multi], present[multi], svtype[multi]


def precedence(lo, hi, present):
    """Per component and caller, whether it holds the union's start / end.

    Two readings, and they answer different questions. *Holds* counts a caller
    that ties for the boundary, which on the 1 kb grid is the common case, and
    says whether the caller is at the edge of the component. *Sole* counts only
    the caller that reaches the boundary alone, which is the coordinate the
    union contributes and an agreement-region policy would discard.
    """
    starts = np.where(present, lo, np.iinfo(np.int64).max)
    ends = np.where(present, hi, np.iinfo(np.int64).min)
    holds_start = (starts == starts.min(axis=1, keepdims=True)) & present
    holds_end = (ends == ends.max(axis=1, keepdims=True)) & present
    sole_start = holds_start & (holds_start.sum(axis=1, keepdims=True) == 1)
    sole_end = holds_end & (holds_end.sum(axis=1, keepdims=True) == 1)
    return holds_start, holds_end, sole_start, sole_end


# --------------------------------------------------------------------------- #
# The table: precedence and span, every coverage
# --------------------------------------------------------------------------- #
rows = []
extents = {}
for coverage in COVERAGES:
    lo, hi, present, svtype = per_caller_extents(coverage)
    extents[coverage] = (lo, hi, present, svtype)
    holds_start, holds_end, sole_start, sole_end = precedence(lo, hi, present)
    union_span = (
        np.where(present, hi, np.iinfo(np.int64).min).max(axis=1)
        - np.where(present, lo, np.iinfo(np.int64).max).min(axis=1)
    )
    for index, caller in enumerate(CALLERS):
        here = present[:, index]
        own_span = hi[here, index] - lo[here, index]
        rows.append({
            "coverage": coverage,
            "caller": LABELS[caller],
            "n_components": int(here.sum()),
            "pct_holds_start": 100.0 * holds_start[here, index].mean(),
            "pct_holds_end": 100.0 * holds_end[here, index].mean(),
            "pct_sole_start": 100.0 * sole_start[here, index].mean(),
            "pct_sole_end": 100.0 * sole_end[here, index].mean(),
            "pct_holds_neither": 100.0 * (~holds_start[here, index] & ~holds_end[here, index]).mean(),
            "median_own_over_union_span": float(np.median(own_span / union_span[here])),
            "median_own_span": float(np.median(own_span)),
        })

table = pd.DataFrame(rows)
TABLES.mkdir(parents=True, exist_ok=True)
table.to_csv(TABLES / "breakpoint_precedence.csv", index=False)

# Pairwise outward extension at the focus coverage: how far beyond the other
# caller each one reaches, positive outward, one row per component and pair.
lo, hi, present, svtype = extents[FOCUS]
pairs = []
for first, second in itertools.combinations(range(len(CALLERS)), 2):
    both = present[:, first] & present[:, second]
    left = lo[both, second] - lo[both, first]   # > 0: `first` starts earlier
    right = hi[both, first] - hi[both, second]  # > 0: `first` ends later
    pairs.append({
        "coverage": FOCUS,
        "caller": LABELS[CALLERS[first]],
        "against": LABELS[CALLERS[second]],
        "n_components": int(both.sum()),
        "pct_identical_both_edges": 100.0 * ((left == 0) & (right == 0)).mean(),
        "pct_extends_left": 100.0 * (left > 0).mean(),
        "pct_extends_right": 100.0 * (right > 0).mean(),
        "median_left_bp": float(np.median(left)),
        "median_right_bp": float(np.median(right)),
        "median_span_ratio": float(np.median((hi[both, first] - lo[both, first])
                                             / (hi[both, second] - lo[both, second]))),
        "pct_strictly_wider": 100.0 * ((hi[both, first] - lo[both, first])
                                       > (hi[both, second] - lo[both, second])).mean(),
    })
    pairs[-1]["_left"], pairs[-1]["_right"] = left, right

pairwise = pd.DataFrame(pairs).drop(columns=["_left", "_right"])
pairwise.to_csv(TABLES / "breakpoint_precedence_pairwise.csv", index=False)

# By variant class at the focus coverage: sole precedence per caller, and the
# share of components in which the caller alone supplies at least one of the
# two union boundaries -- the components whose consensus coordinates would
# change if that caller's calls were trimmed back to the others'.
_, _, sole_start, sole_end = precedence(lo, hi, present)
by_class = pd.DataFrame([
    {
        "coverage": FOCUS,
        "class": cls,
        "caller": LABELS[caller],
        "n_components_in_class": int((svtype == cls).sum()),
        "n_components": int(((svtype == cls) & present[:, index]).sum()),
        "pct_sole_start": 100.0 * sole_start[(svtype == cls) & present[:, index], index].mean(),
        "pct_sole_end": 100.0 * sole_end[(svtype == cls) & present[:, index], index].mean(),
        "pct_class_with_sole_boundary": 100.0 * (sole_start[:, index] | sole_end[:, index])[svtype == cls].mean(),
    }
    for cls in ("DEL", "DUP")
    for index, caller in enumerate(CALLERS)
])
by_class.to_csv(TABLES / "breakpoint_precedence_by_class.csv", index=False)


# --------------------------------------------------------------------------- #
# The figure
# --------------------------------------------------------------------------- #
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(7.09, 2.7), width_ratios=[1.25, 1.0])
fig.subplots_adjust(left=0.115, right=0.985, top=0.88, bottom=0.235, wspace=0.34)


def panel(ax, letter: str, x: float) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=3.2)
    ax.text(x, 1.04, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")


# (A) Pairwise outward extension. Each pair contributes two boxes, one per edge,
# and the sign says which caller reaches further: above zero the named caller
# extends beyond the one it is compared with. Whiskers at the 5th and 95th
# percentiles, since the tails run to hundreds of kilobases on a few components.
positions, boxes, tick_labels, tick_positions = [], [], [], []
for index, row in enumerate(pairs):
    for offset, edge in ((-0.17, "_left"), (0.17, "_right")):
        positions.append(index + offset)
        boxes.append(row[edge])
    tick_positions.append(index)
    tick_labels.append(f"{row['caller']}\nvs {row['against']}\n"
                       f"({row['pct_identical_both_edges']:.0f}% identical)")

drawn = ax_a.boxplot(boxes, positions=positions, widths=0.28, showfliers=False,
                     whis=(5, 95), patch_artist=True)
for index, patch in enumerate(drawn["boxes"]):
    patch.set_facecolor(COLORS[pairs[index // 2]["caller"]])
    patch.set_alpha(0.35 if index % 2 else 0.85)
    patch.set_linewidth(0.6)
for key in ("medians", "whiskers", "caps"):
    for artist in drawn[key]:
        artist.set_color("#222222")
        artist.set_linewidth(0.7)
ax_a.axhline(0, color="#9E9E9E", linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
ax_a.set_xticks(tick_positions, tick_labels)
ax_a.set_ylim(-2400, 3100)
ax_a.set_ylabel("Outward extension of the\nnamed caller (bp)", labelpad=2)
ax_a.legend(handles=[
    Patch(facecolor="#767676", alpha=0.85, label="Start (5' edge)"),
    Patch(facecolor="#767676", alpha=0.35, label="End (3' edge)"),
], frameon=False, loc="upper center", ncol=2, handlelength=1.1, borderpad=0,
    columnspacing=1.0, handletextpad=0.5)
panel(ax_a, "A", x=-0.15)

# (B) Boundary precedence by variant class: of the components a caller appears
# in, the share in which it *alone* reaches the union start or the union end.
# The sole reading rather than the tie-inclusive one, because a boundary two
# callers agree on is not a coordinate the union had to go outside them to get.
_, _, holds_start, holds_end = precedence(lo, hi, present)
width = 0.26
offsets = {"DEL": -width / 2, "DUP": width / 2}
# Two class groups per caller, each holding the two edges. Class is labelled on
# a minor tick row so the caller names can sit under the group they label.
centres = np.arange(len(CALLERS), dtype=float)
minor_positions, minor_labels = [], []
for cls, offset in (("DEL", -0.19), ("DUP", 0.19)):
    rows_cls = svtype == cls
    for edge, (holds, alpha) in enumerate(((holds_start, 0.9), (holds_end, 0.4))):
        shares = [
            100.0 * holds[rows_cls & present[:, index], index].mean()
            for index in range(len(CALLERS))
        ]
        ax_b.bar(centres + offset + (edge - 0.5) * 0.155, shares, width=0.15, alpha=alpha,
                 color=[COLORS[LABELS[c]] for c in CALLERS],
                 edgecolor="white", linewidth=0.4, zorder=2)
    minor_positions.extend(centres + offset)
    minor_labels.extend([cls] * len(CALLERS))
ax_b.set_xticks(centres, [LABELS[c] for c in CALLERS])
ax_b.set_xticks(minor_positions, minor_labels, minor=True)
ax_b.tick_params(axis="x", which="major", length=0, pad=13)
ax_b.tick_params(axis="x", which="minor", length=0, labelsize=6, pad=2)
ax_b.set_xlim(-0.55, len(CALLERS) - 0.45)
ax_b.set_ylim(0, 70)
ax_b.set_ylabel("Components whose union boundary\nthe caller alone reaches (%)", labelpad=2)
ax_b.legend(handles=[
    Patch(facecolor="#767676", alpha=0.9, label="Start (5' edge)"),
    Patch(facecolor="#767676", alpha=0.4, label="End (3' edge)"),
], frameon=False, loc="upper center", ncol=2, handlelength=1.1, borderpad=0,
    columnspacing=1.0, handletextpad=0.5)
panel(ax_b, "B", x=-0.24)

DEST.mkdir(parents=True, exist_ok=True)
for suffix, dpi in ((".png", 600), (".pdf", None)):
    fig.savefig(DEST / f"breakpoint_precedence{suffix}", dpi=dpi)
plt.close(fig)

print(table.round(1).to_string(index=False))
print()
print(pairwise.round(2).to_string(index=False))
print()
print(by_class.round(1).to_string(index=False))
print(f"\nwrote {TABLES / 'breakpoint_precedence.csv'}, "
      f"{TABLES / 'breakpoint_precedence_pairwise.csv'} and "
      f"{DEST / 'breakpoint_precedence.png'}")
