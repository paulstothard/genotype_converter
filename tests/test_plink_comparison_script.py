from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_compare_module():
    path = Path("validation/plink_comparison/scripts/compare_plink_bim.py")
    spec = importlib.util.spec_from_file_location("compare_plink_bim", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["compare_plink_bim"] = module
    spec.loader.exec_module(module)
    return module


def test_compare_plink_bim_classifies_allele_relationships(tmp_path):
    module = _load_compare_module()
    expected = tmp_path / "expected.bim"
    actual = tmp_path / "actual.bim"
    expected.write_text(
        "\n".join(
            [
                "1 exact 0 10 A C",
                "1 swapped 0 20 A G",
                "1 complement 0 30 A C",
                "1 mismatch 0 40 A C",
                "1 missing 0 50 A C",
            ]
        )
        + "\n"
    )
    actual.write_text(
        "\n".join(
            [
                "1 exact 0 10 A C",
                "1 swapped 0 20 G A",
                "1 complement 0 30 T G",
                "1 mismatch 0 41 A G",
                "1 extra 0 60 A C",
            ]
        )
        + "\n"
    )

    rows = {row["marker_name"]: row for row in module.compare(expected, actual)}

    assert rows["exact"]["status"] == "exact"
    assert rows["swapped"]["status"] == "swapped"
    assert rows["complement"]["status"] == "complement"
    assert rows["mismatch"]["status"] == "mismatch"
    assert rows["mismatch"]["position_status"] == "different"
    assert rows["missing"]["status"] == "missing_from_genotype_converter_plus"
    assert rows["extra"]["status"] == "extra_in_genotype_converter_plus"
