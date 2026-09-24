"""The exclusion policy: drop a call whole past `max_excluded_fraction`, else keep
it whole, trimming only its ends out of the mask. Never split a call.

Importing `consensuscnv.parsing` loads cyvcf2 and pandas (~500 ms), so this is the
one slow module in the suite.
"""

import pytest

from consensuscnv.parsing.parser_utils import ExclusionMask
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
