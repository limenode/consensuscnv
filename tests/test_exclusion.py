"""The exclusion policy: drop a call whole past `max_excluded_fraction`, else keep
it whole, trimming only its ends out of the mask. Never split a call.

The hand-written cases pin the semantics; the randomized ones check the interval
arithmetic against `reference_apply`, which restates the policy one base at a
time, on masks with overlapping and touching rows and calls whose ends sit on and
around region edges.

Importing `consensuscnv.parsing` loads cyvcf2 and pandas (~500 ms), so this is the
one slow module in the suite.
"""

import random

import pytest

from consensuscnv.parsing.parser_utils import ExclusionMask
from consensuscnv.parsing.vcf_parser import _process_single_vcf_to_df
from consensuscnv.utils import PipelineConfig

# chr1: [1000, 2000) and [3000, 3500) written as two touching halves, so loading
# must merge them for the fraction to count 500 bp once; [10_000, 10_100) is an
# island well inside the long calls below.
MASK_ROWS = (
    "# header lines are skipped\n"
    "chr1\t1000\t2000\n"
    "chr1\t3000\t3250\n"
    "1\t3250\t3500\n"
    "chr1\t10000\t10100\n"
)


@pytest.fixture
def mask(tmp_path) -> ExclusionMask:
    path = tmp_path / "mask.bed"
    path.write_text(MASK_ROWS)
    return ExclusionMask.load(path)


def test_load_merges_touching_regions_and_adds_chr_prefix(mask):
    assert mask.regions == {"chr1": [(1000, 2000), (3000, 3500), (10_000, 10_100)]}


def test_fraction_sums_every_overlapping_region(mask):
    # [1500, 3500): 500 bp from the first region plus 500 from the second, of 2000.
    assert mask.excluded_fraction("chr1", 1500, 3500) == 0.5
    assert not mask.is_excluded("chr1", 1500, 3500, 0.5)  # strictly more than
    assert mask.is_excluded("chr1", 1500, 3500, 0.49)


def test_interior_island_is_kept_inside_an_unsplit_call(mask):
    assert mask.apply("chr1", 5000, 20_000, 0.5, trim_ends=True) == (5000, 20_000)


def test_ends_inside_the_mask_are_trimmed_to_its_edges(mask):
    # Start inside [1000, 2000), end inside [3000, 3500): 1000 of 2500 bp masked.
    assert mask.apply("chr1", 1500, 3200, 0.5, trim_ends=True) == (2000, 3000)
    assert mask.apply("chr1", 1500, 3200, 0.5, trim_ends=False) == (1500, 3200)


def test_trimming_leaves_interior_regions_alone(mask):
    # Start trimmed out of [1000, 2000); [3000, 3500) and the island stay inside.
    assert mask.apply("chr1", 1900, 12_000, 0.5, trim_ends=True) == (2000, 12_000)


def test_call_past_the_fraction_is_dropped_whole(mask):
    assert mask.apply("chr1", 900, 2100, 0.5, trim_ends=True) is None


def test_call_inside_one_region_is_dropped_even_at_fraction_one(mask):
    # Kept by the fraction (1.0 is not more than 1.0), then empty once trimmed.
    assert mask.apply("chr1", 1200, 1800, 1.0, trim_ends=False) == (1200, 1800)
    assert mask.apply("chr1", 1200, 1800, 1.0, trim_ends=True) is None


def test_unmasked_chromosome_is_untouched(mask):
    assert mask.apply("chr2", 1000, 2000, 0.0, trim_ends=True) == (1000, 2000)


def _config(tmp_path, **raw) -> PipelineConfig:
    genome = tmp_path / "genome.txt"
    genome.write_text("chr1\t1000000\n")
    return PipelineConfig.from_raw(
        {
            "experimental": {"TS": {"caller": "/nowhere/{id}.vcf"}},
            "output_dir": str(tmp_path / "out"),
            "genome_file": str(genome),
            **raw,
        }
    )


def test_config_defaults_keep_calls_whole_and_trim_their_ends(tmp_path):
    config = _config(tmp_path)
    assert config.max_excluded_fraction == 0.5
    assert config.trim_excluded_ends is True


def test_config_rejects_a_non_boolean_trim_flag(tmp_path):
    with pytest.raises(ValueError, match="trim_excluded_ends"):
        _config(tmp_path, trim_excluded_ends="no")


# --------------------------------------------------------------------------- #
# Randomized: the package against a base-by-base restatement of the policy.
# --------------------------------------------------------------------------- #

CHROM_LENGTH = 3000
SEEDS = range(8)


def reference_apply(masked, start, end, max_excluded_fraction, trim_ends):
    """`ExclusionMask.apply`, one base at a time: `masked` is a set of positions."""
    covered = sum(base in masked for base in range(start, end))
    if covered > max_excluded_fraction * (end - start):
        return None
    if trim_ends:
        while start < end and start in masked:
            start += 1
        while end > start and end - 1 in masked:
            end -= 1
        if start >= end:
            return None
    return start, end


def random_case(seed, n_rows=25, n_calls=400):
    """Unsorted, overlapping and touching mask rows, and calls whose ends are
    drawn mostly from region edges and their neighbours."""
    rng = random.Random(seed)
    rows = []
    for _ in range(n_rows):
        start = rng.randrange(CHROM_LENGTH - 1)
        rows.append((start, min(CHROM_LENGTH, start + rng.randint(1, 150))))
    rows += [(end, end + 20) for _, end in rng.sample(rows, 3) if end + 20 <= CHROM_LENGTH]

    edges = sorted({p + d for row in rows for p in row for d in (-1, 0, 1)} - {-1})
    calls = []
    while len(calls) < n_calls:
        pick = [rng.choice(edges) if rng.random() < 0.7 else rng.randrange(CHROM_LENGTH + 1)
                for _ in range(2)]
        start, end = min(pick), max(pick)
        if start < end <= CHROM_LENGTH:
            calls.append((start, end))

    masked = {base for start, end in rows for base in range(start, end)}
    return rows, calls, masked


def load_mask(tmp_path, rows) -> ExclusionMask:
    path = tmp_path / "mask.bed"
    path.write_text("".join(f"chr1\t{start}\t{end}\n" for start, end in rows))
    return ExclusionMask.load(path)


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("fraction", [0.0, 0.25, 0.5, 1.0])
@pytest.mark.parametrize("trim", [True, False])
def test_apply_matches_the_base_by_base_reference(tmp_path, seed, fraction, trim):
    rows, calls, masked = random_case(seed)
    mask = load_mask(tmp_path, rows)
    for start, end in calls:
        covered = sum(base in masked for base in range(start, end))
        assert mask.overlap_bp("chr1", start, end) == covered, (start, end)
        assert mask.apply("chr1", start, end, fraction, trim) == reference_apply(
            masked, start, end, fraction, trim
        ), (start, end)


@pytest.mark.parametrize("seed", SEEDS)
def test_no_kept_call_is_over_half_masked_or_ends_inside_the_mask(tmp_path, seed):
    rows, calls, masked = random_case(seed)
    mask = load_mask(tmp_path, rows)
    n_kept = n_holding_islands = 0
    for start, end in calls:
        kept = mask.apply("chr1", start, end, 0.5, trim_ends=True)
        if kept is None:
            continue
        n_kept += 1
        # Measured on the call as reported, before trimming.
        assert sum(base in masked for base in range(start, end)) <= 0.5 * (end - start)
        new_start, new_end = kept
        assert start <= new_start < new_end <= end
        # Neither the first nor the last base of a kept call is masked ...
        assert new_start not in masked and new_end - 1 not in masked, (start, end, kept)
        # ... so every mask region it touches lies wholly inside it: never split.
        for region_start, region_end in mask.overlapping("chr1", new_start, new_end):
            assert new_start < region_start and region_end < new_end
            n_holding_islands += 1
    # Guard against a vacuous pass: the cases must exercise both properties.
    assert n_kept > 0 and n_holding_islands > 0


@pytest.mark.parametrize("seed", SEEDS)
def test_vcf_parser_output_and_statistics_match_the_reference(tmp_path, seed):
    """End to end through `_process_single_vcf_to_df`: the BED rows written and
    every exclusion counter, against the reference."""
    rows, calls, masked = random_case(seed)
    mask = load_mask(tmp_path, rows)
    calls = sorted(calls)

    vcf = tmp_path / "calls.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        f"##contig=<ID=chr1,length={CHROM_LENGTH}>\n"
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="type">\n'
        '##INFO=<ID=END,Number=1,Type=Integer,Description="end">\n'
        '##ALT=<ID=DEL,Description="deletion">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(
            f"chr1\t{start + 1}\t.\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END={end}\n"
            for start, end in calls
        )
    )
    df, stats = _process_single_vcf_to_df(
        vcf, mask, ("chr1",), max_excluded_fraction=0.5, trim_excluded_ends=True
    )

    expected, dropped, trimmed = [], [], []
    for start, end in calls:
        kept = reference_apply(masked, start, end, 0.5, True)
        if kept is None:
            dropped.append((start, end))
        else:
            expected.append(kept)
            if kept != (start, end):
                trimmed.append((end - start) - (kept[1] - kept[0]))

    assert list(zip(df["start"], df["end"])) == expected
    assert stats["total_call_count"] == len(calls)
    assert stats["calls_removed_excluded"] == stats["calls_del_removed_excluded"] == len(dropped)
    assert stats["bases_removed_excluded"] == sum(end - start for start, end in dropped)
    assert stats["bases_masked_excluded"] == sum(
        sum(base in masked for base in range(start, end)) for start, end in dropped
    )
    assert stats["calls_trimmed_excluded"] == len(trimmed)
    assert stats["bases_trimmed_excluded"] == sum(trimmed)
    assert dropped and trimmed  # both paths exercised
