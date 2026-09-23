"""How much the sequencing-vs-array comparison leans on the shared-library benchmark source.

Where used
----------
Discussion -> the limitations paragraph on benchmark independence: the
1000 Genomes high-coverage SV call set was called from the same 30x alignments
the sequencing call sets were derived from, so its errors and blind spots are
shared with the query and not with the SNP array.
    Supplementary Note S4                        every number in the note
    Supplementary Table Benchmark Independence   precision and recall against
                                                 each benchmark, by class

The concern is concrete. A benchmark source called from the same reads can
confirm an artifact both pipelines produce from those reads, which would score a
sequencing call as a true positive for a reason the array has no access to. Two
measurements bound how much that could matter:

  * *reliance*: of the benchmark intervals each call set recovers against the
    full merged benchmark, the share supported by the 1000 Genomes set alone.
    Recoveries of this kind are the only ones the shared library could have
    produced, so the share is an upper bound on its contribution to recall.
  * *re-scoring*: every call set scored again against a benchmark merged from
    the two long-read sources only (HGSVC3 and ONT Vienna), neither of which
    shares reads with the query. If the sequencing call sets' advantage over the
    array were an artifact of the shared library it would shrink or vanish here.

Deletions are the comparison to read. HGSVC3 carries no duplications and ONT
Vienna's come from SVAN-annotated insertions, so a long-read-only benchmark
holds far fewer duplications than the full one and its duplication metrics
describe a different truth set rather than the same one without 1000G.

Adopted parameters throughout: padding 0, 1 kb floor on both sides, 2-of-3
consensus at 0.5, classification at 0.5.

    pixi run python manuscript/scripts/benchmark_independence.py
"""

import glob
from pathlib import Path

import pandas as pd

from consensuscnv.callsets import collect_callsets, merge_components
from consensuscnv.callsets.registry import SOURCES, SVTYPES, seed_chromosomes
from consensuscnv.classification.classify import classify
from consensuscnv.classification.intervals import IntervalSet
from consensuscnv.classification.pairs import build_candidates
from consensuscnv.utils import read_genome_file

ROOT = Path("/lab01/Projects/Lionel_Projects/blendedCNV_pipeline")

# Chromosome ids have to order the genome rather than the order the BEDs
# happen to be read in, so the registry is seeded before any CallSet is built.
seed_chromosomes(read_genome_file(ROOT / "src" / "consensuscnv" / "templates" / "genome_primary_hg38.txt"))
TABLES = ROOT / "results" / "manuscript"

CALLERS = ("cnvpytor", "delly", "gatk")
COVERAGES = ("30x", "6x", "4x", "2x")
ARRAY = "SNP array"
ORDER = (*COVERAGES, ARRAY)
SAMPLES = (
    "HG00096", "HG00171", "HG00268", "HG00513", "HG00731", "HG01596", "HG01890",
    "NA18989", "NA19129", "NA19238", "NA19331", "NA19347", "NA20847",
)
BENCHMARKS = {
    "full": ("1000G", "HGSVC3", "ont_vienna"),
    "long-read only": ("HGSVC3", "ont_vienna"),
}
SHARED_SOURCE = "1000g"  # the label in BED column 5, not the directory name

SIZE_FLOOR = 1_000
BENCHMARK_PADDING = 0
CONSENSUS_THRESHOLD = 0.5
CLASSIFY_THRESHOLD = 0.5
CONSENSUS_LEVEL = 2


def bed_paths(root: Path, subdirs: tuple[str, ...]) -> list[str]:
    """Every per-sample BED under the named subdirectories, never a wildcard."""
    return [bed for sub in subdirs for bed in sorted(glob.glob(str(root / sub / "*.bed")))]


def floored(intervals: IntervalSet) -> IntervalSet:
    return intervals.filter_by_size(min_size=SIZE_FLOOR)


def benchmark(sources: tuple[str, ...]) -> IntervalSet:
    return floored(IntervalSet.from_merged(merge_components(
        collect_callsets(bed_paths(ROOT / "out" / "benchmark", sources)),
        max_padding=BENCHMARK_PADDING,
    )))


def consensus_at(coverage: str) -> IntervalSet:
    """The 2-of-3 set, off one merge selected on n_sources (== min_sources=2)."""
    merged = IntervalSet.from_merged(merge_components(
        collect_callsets(bed_paths(ROOT / "out" / f"{coverage}_Coverage", CALLERS)),
        min_reciprocal_overlap=CONSENSUS_THRESHOLD,
    ))
    return floored(merged.select(merged.n_sources >= CONSENSUS_LEVEL))


truths = {name: benchmark(sources) for name, sources in BENCHMARKS.items()}
queries = {
    **{coverage: consensus_at(coverage) for coverage in COVERAGES},
    ARRAY: floored(IntervalSet.from_bed(ROOT / "out" / "SNP_Array" / "bed" / f"{s}.bed" for s in SAMPLES)),
}
CLASS_IDS = {"DEL": SVTYPES.intern("DEL"), "DUP": SVTYPES.intern("DUP")}
shared_bit = 1 << SOURCES.intern(SHARED_SOURCE)


def restricted(intervals: IntervalSet, cls: str | None) -> IntervalSet:
    return intervals if cls is None else intervals.select(intervals.svtype_idx == CLASS_IDS[cls])


rows = []
for benchmark_name, truth_all in truths.items():
    for cls in (None, "DEL", "DUP"):
        truth = restricted(truth_all, cls)
        for name in ORDER:
            query = restricted(queries[name], cls)
            result = classify(build_candidates(query, truth),
                              min_reciprocal_overlap=CLASSIFY_THRESHOLD, validate=False)
            summary = result.summary()
            found = result.truth_matched
            row = {
                "benchmark": benchmark_name,
                "class": cls or "all",
                "call set": name,
                "n_query": summary.n_query,
                "n_truth": summary.n_truth,
                "n_true_positive": summary.n_true_positive,
                "n_truth_found": summary.n_truth_found,
                "precision": summary.precision,
                "recall": summary.recall,
            }
            if benchmark_name == "full":
                # Recoveries whose benchmark interval rests on the shared-library
                # source alone: the only ones that source could have produced.
                only_shared = truth.source_bits == shared_bit
                row["n_found_shared_only"] = int((found & only_shared).sum())
                row["pct_found_shared_only"] = 100.0 * (found & only_shared).sum() / max(found.sum(), 1)
                row["pct_truth_shared_only"] = 100.0 * only_shared.mean()
            rows.append(row)

table = pd.DataFrame(rows)
TABLES.mkdir(parents=True, exist_ok=True)
table.to_csv(TABLES / "benchmark_independence.csv", index=False)

pd.set_option("display.width", 200)
print(table.round(3).to_string(index=False))
print(f"\nwrote {TABLES / 'benchmark_independence.csv'}")
