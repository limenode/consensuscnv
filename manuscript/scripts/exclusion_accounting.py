"""What the exclusion mask does to each caller's calls and to the consensus sets.

Where used
----------
Supplementary Tables -> "Supplementary Table Exclusion Accounting" (per caller)
and "Supplementary Table Exclusion Consensus" (per consensus level), cited from
Methods -> "Parsing" as the comparison of mask-overlap stringencies.

Two policies, both applied through the package's own `ExclusionMask.apply`:

    adopted   max_excluded_fraction 0.5, trim_excluded_ends true  (internal.config.yaml)
    strict    max_excluded_fraction 0.01, trim_excluded_ends false (the earlier default)

Per caller, every raw call (parsed with no mask) falls in exactly one class:
dropped (masked beyond the threshold); kept with an end trimmed out of the mask;
kept unchanged but sharing sequence with the mask ("overlap": under the adopted
policy only whole regions inside the call, under the strict one also an end
reaching up to 1% into it); or kept with no overlap at all. The
kept calls are checked call for call against the package's own masked parse.

Per consensus level, the unmasked and masked merges are linked through the calls
each masked component came from, and N(no mask) - N(mask) is split exactly into:
vanished (every member dropped), lost a caller (fell below the level), below the
floor (trimmed under `min_size`), gained, and a net for components that split or
merged. Consensus uses the config's `consensus:` block (adopted 0.5 reciprocal
overlap, 1 kb floor), as `consensuscnv call` does.

Nothing is written under `out/`; VCFs are parsed in memory. As a check that the
script reproduces the paper's parse, the strict policy's kept calls are compared
with the BEDs currently under `out/`, which were parsed at 0.01 without trimming.

Outputs, in `results/manuscript/`:
    exclusion_accounting_calls.csv
    exclusion_accounting_consensus.csv
    supp_exclusion_accounting_tables.typ   both tables, ready to paste

    pixi run python manuscript/scripts/exclusion_accounting.py
"""

import glob
from collections import Counter
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from consensuscnv.callsets import Call, collect_callsets, merge_components, read_bed_calls
from consensuscnv.callsets.registry import seed_chromosomes
from consensuscnv.parsing.parser_utils import ExclusionMask
from consensuscnv.parsing.vcf_parser import (
    _process_single_vcf_to_df,
    get_experimental_sets_from_config,
)
from consensuscnv.utils import build_config, load_sample_list

ROOT = Path("/lab01/Projects/Lionel_Projects/blendedCNV_pipeline")
CONFIG = ROOT / "internal.config.yaml"
DEST = ROOT / "results" / "manuscript"

config = build_config(CONFIG)
POLICIES = {
    "adopted": (config.max_excluded_fraction, config.trim_excluded_ends),
    "strict": (0.01, False),
}
assert POLICIES["adopted"] == (0.5, True), "internal.config.yaml no longer pins the adopted policy"

CALLER_NAMES = {"cnvpytor": "CNVpytor", "delly": "Delly", "gatk": "GATK-gCNV"}
MASK = ExclusionMask.load(config.excluded_regions_file)
NO_MASK = ExclusionMask({})
OVERLAP = config.consensus.reciprocal_overlaps
MIN_SIZE = config.consensus.min_size
seed_chromosomes(config.chromosomes)


def parse(path, mask, fraction, trim) -> list[tuple[str, int, int, str]]:
    df, _ = _process_single_vcf_to_df(
        path, mask, config.chromosomes,
        max_excluded_fraction=fraction, trim_excluded_ends=trim,
    )
    return list(df.itertuples(index=False, name=None))


def per_label(merged):
    """n_sources and span per component label, from an unfiltered merge."""
    n = int(merged.labels.max()) + 1
    sources, span = np.zeros(n, np.int64), np.zeros(n, np.int64)
    sources[merged.component_id] = merged.n_sources
    span[merged.component_id] = merged.ends - merged.starts
    return sources, span


def reconcile(calls0, calls1, origin, overlap, n_tools):
    """Split N(no mask) - N(mask) into exact parts, per level, at one overlap."""
    cs0 = collect_callsets([calls0], chromosome_order=config.chromosomes)
    cs1 = collect_callsets([calls1], chromosome_order=config.chromosomes)
    index0 = {c: i for i, c in enumerate(cs0.calls)}  # identical calls share a component
    # Two raw calls can trim to the same masked call: give each copy its own origin.
    pool = {k: list(v) for k, v in origin.items()}
    node0 = np.fromiter((index0[pool[c].pop()] for c in cs1.calls), np.int64, len(cs1.calls))

    m0 = merge_components(cs0, min_reciprocal_overlap=overlap)
    m1 = merge_components(cs1, min_reciprocal_overlap=overlap)
    src0, span0 = per_label(m0)
    src1, span1 = per_label(m1)
    n0, n1 = len(src0), len(src1)

    # Families: connected pieces of the links between unmasked and masked components.
    pairs = np.unique(np.stack([m0.labels[node0], m1.labels]), axis=1)
    graph = coo_matrix(
        (np.ones(pairs.shape[1]), (pairs[0], n0 + pairs[1])), shape=(n0 + n1, n0 + n1)
    )
    _, family = connected_components(graph, directed=False)
    fam0, fam1 = family[:n0], family[n0:]
    size0 = np.bincount(fam0, minlength=family.max() + 1)
    size1 = np.bincount(fam1, minlength=family.max() + 1)
    one0 = (size0[fam0] == 1) & (size1[fam0] == 1)
    one1 = (size0[fam1] == 1) & (size1[fam1] == 1)
    partner = np.zeros(family.max() + 1, np.int64)
    partner[fam0[one0]] = np.flatnonzero(one0)
    p1 = partner[fam1]  # the unmasked partner of each one-to-one masked component

    rows = []
    for level in range(1, n_tools + 1):
        q0 = (src0 >= level) & (span0 >= MIN_SIZE)
        q1 = (src1 >= level) & (span1 >= MIN_SIZE)
        down = one1 & q0[p1] & ~q1
        row = {
            "level": level,
            "no_mask": int(q0.sum()),
            "vanished": int(q0[size1[fam0] == 0].sum()),
            "lost_caller": int((down & (src1 < level)).sum()),
            "below_floor": int((down & (src1 >= level)).sum()),
            "gained": int((one1 & ~q0[p1] & q1).sum()),
            "split_merge_net": int((q0 & (size1[fam0] > 0) & ~one0).sum() - (q1 & ~one1).sum()),
            "masked": int(q1.sum()),
        }
        assert row["no_mask"] - row["masked"] == (
            row["vanished"] + row["lost_caller"] + row["below_floor"]
            - row["gained"] + row["split_merge_net"]
        ), row
        rows.append(row)
    return rows


samples = load_sample_list(config.sample_list_file)
call_rows, consensus_rows = [], []

for call_set, tools in get_experimental_sets_from_config(config, samples).items():
    coverage = call_set.split()[0]
    # As `process_vcfs_to_beds`: only samples every tool has.
    common = set.intersection(*(set(m) for m in tools.values()))
    raw = {tool: [] for tool in tools}
    for tool, sample_map in tools.items():
        for sample_id in sorted(common):
            raw[tool] += [
                Call(chrom, start, end, svtype, tool, sample_id)
                for chrom, start, end, svtype in parse(sample_map[sample_id], NO_MASK, 1.0, False)
            ]
    calls0 = [c for tool in tools for c in raw[tool]]

    for policy, (fraction, trim) in POLICIES.items():
        calls1, origin = [], {}
        for tool, sample_map in tools.items():
            n = Counter()
            for c in raw[tool]:
                kept = MASK.apply(c.chrom, c.start, c.end, fraction, trim)
                if kept is None:
                    n["dropped"] += 1
                    continue
                k = replace(c, start=kept[0], end=kept[1])
                if kept != (c.start, c.end):
                    n["trimmed"] += 1
                elif MASK.overlap_bp(c.chrom, c.start, c.end):
                    n["island"] += 1
                else:
                    n["untouched"] += 1
                calls1.append(k)
                origin.setdefault(k, []).append(c)

            # The prediction must equal the package's own masked parse, call for call.
            parsed = Counter(
                Call(chrom, start, end, svtype, tool, sample_id)
                for sample_id in sorted(common)
                for chrom, start, end, svtype in parse(sample_map[sample_id], MASK, fraction, trim)
            )
            predicted = Counter(k for k in calls1 if k.source == tool)
            assert predicted == parsed, f"{call_set} {tool} {policy}: prediction != parser"
            if policy == "strict":
                on_disk = Counter(
                    replace(c, source=tool)
                    for p in glob.glob(str(config.layout.bed_tool_dir(call_set, tool) / "*.bed"))
                    for c in read_bed_calls(p)
                )
                print(f"  {call_set} {tool}: strict policy "
                      f"{'matches' if on_disk == predicted else 'DIFFERS FROM'} out/")

            call_rows.append({
                "coverage": coverage, "caller": tool, "policy": policy,
                "raw": len(raw[tool]), **{k: n[k] for k in ("dropped", "trimmed", "island", "untouched")},
                "kept": sum(predicted.values()),
            })

        for overlap in OVERLAP:
            for row in reconcile(calls0, calls1, origin, overlap, len(tools)):
                consensus_rows.append(
                    {"coverage": coverage, "policy": policy, "overlap": overlap, **row}
                )

calls = pd.DataFrame(call_rows)
consensus = pd.DataFrame(consensus_rows)
assert (calls.raw == calls.dropped + calls.trimmed + calls.island + calls.untouched).all()
DEST.mkdir(parents=True, exist_ok=True)
calls.to_csv(DEST / "exclusion_accounting_calls.csv", index=False)
consensus.to_csv(DEST / "exclusion_accounting_consensus.csv", index=False)

pd.set_option("display.width", 200)
print(calls.to_string(index=False))
print(consensus.to_string(index=False))


# --------------------------------------------------------------------------- #
# Typst, in the style of the other Supplementary Tables.
# --------------------------------------------------------------------------- #

def num(value) -> str:
    return f"{value:,}".replace("-", "−")


def block(columns, align, header, body) -> list[str]:
    return [
        "#block(width: 100%)[",
        "#set text(hyphenate: false, size: 9.5pt)",
        "#show table.cell.where(y: 1): strong",
        "#table(",
        f"  columns: {columns},",
        f"  align: {align},",
        "  stroke: none,",
        "  table.hline(stroke: 0.6pt),",
        *header,
        "  table.hline(stroke: 0.6pt),",
        *body,
        "  table.hline(stroke: 0.6pt),",
        ")",
        "]",
    ]


wide = calls.pivot_table(index=["coverage", "caller"], columns="policy", values=[
    "raw", "dropped", "trimmed", "island", "kept"]).astype(int)
body, coverages = [], list(dict.fromkeys(calls.coverage))
for i, coverage in enumerate(coverages):
    if i:
        body.append("  table.hline(stroke: 0.3pt),")
    for j, caller in enumerate(CALLER_NAMES):
        r = wide.loc[(coverage, caller)]
        cells = [coverage if j == 0 else "", CALLER_NAMES[caller], num(r[("raw", "adopted")]),
                 num(r[("dropped", "strict")]), num(r[("island", "strict")]), num(r[("kept", "strict")]),
                 num(r[("dropped", "adopted")]), num(r[("trimmed", "adopted")]),
                 num(r[("island", "adopted")]), num(r[("kept", "adopted")])]
        body.append("  " + ", ".join(f"[{c}]" for c in cells) + ",")
calls_typ = block(
    "(auto, 1fr, auto, auto, auto, auto, auto, auto, auto, auto)",
    "(left, left, right, right, right, right, right, right, right, right)",
    [
        "  table.header(",
        "    table.cell(rowspan: 2, align: left + bottom)[Coverage],",
        "    table.cell(rowspan: 2, align: left + bottom)[Caller],",
        "    table.cell(rowspan: 2, align: right + bottom)[Raw Calls],",
        "    table.cell(colspan: 3, align: center)[Strict (1%)],",
        "    table.cell(colspan: 4, align: center)[Adopted (50%, trimmed)],",
        "    table.hline(start: 3, end: 6, stroke: 0.3pt),",
        "    table.hline(start: 6, end: 10, stroke: 0.3pt),",
        "    [Dropped], [Overlap], [Kept], [Dropped], [Trimmed], [Overlap], [Kept],",
        "  ),",
    ],
    body,
)

body = []
for i, coverage in enumerate(coverages):
    if i:
        body.append("  table.hline(stroke: 0.3pt),")
    sub = consensus[consensus.coverage == coverage].set_index(["policy", "level"])
    for level in (1, 2, 3):
        a, s = sub.loc[("adopted", level)], sub.loc[("strict", level)]
        cells = [coverage if level == 1 else "", f"{level}/3", num(a.no_mask),
                 num(a.vanished), num(a.lost_caller), num(a.below_floor), num(a.gained),
                 num(a.split_merge_net), num(a.masked), num(s.masked)]
        body.append("  " + ", ".join(f"[{c}]" for c in cells) + ",")
consensus_typ = block(
    "(auto, auto, auto, auto, auto, auto, auto, auto, auto, auto)",
    "(left, left, right, right, right, right, right, right, right, right)",
    [
        "  table.header(",
        "    table.cell(rowspan: 2, align: left + bottom)[Coverage],",
        "    table.cell(rowspan: 2, align: left + bottom)[Level],",
        "    table.cell(rowspan: 2, align: right + bottom)[No Mask],",
        "    table.cell(colspan: 6, align: center)[Adopted (50%, trimmed)],",
        "    table.cell(rowspan: 2, align: right + bottom)[Strict (1%)],",
        "    table.hline(start: 3, end: 9, stroke: 0.3pt),",
        "    [Vanished], [Lost Caller], [Below Floor], [Gained], [Split/Merge], [Masked],",
        "  ),",
    ],
    body,
)

(DEST / "supp_exclusion_accounting_tables.typ").write_text(
    "\n".join(calls_typ + [""] + consensus_typ) + "\n"
)
print(f"\nwrote {DEST / 'supp_exclusion_accounting_tables.typ'}")
