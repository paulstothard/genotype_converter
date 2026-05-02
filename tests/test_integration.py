from __future__ import annotations

import csv
from pathlib import Path


def _read_data_rows(csv_path: Path) -> list[dict]:
    """Read CSV skipping comment lines (starting with #)."""
    with open(csv_path) as f:
        lines = [l for l in f if not l.startswith("#")]
    return list(csv.DictReader(lines))


def test_position_matches_expected(pipeline_output):
    rows = _read_data_rows(pipeline_output / "manifest.reference.position.csv")
    by_name = {r["marker_name"]: r for r in rows}

    expected = {
        # SNPs — values validated against original pipeline sample_output
        "SNP1": ("1", "300", "A",     "G"),
        "SNP2": ("1", "700", "T",     "G"),
        "SNP3": ("1", "900", "G",     "T"),
        "SNP4": ("1", "200", "A",     "G"),
        "SNP5": ("1", "800", "T",     "G"),
        # Indels — anchor-based VCF positions computed from artificial reference
        # INDEL1: insertion [-/ACGT], anchor at ref pos 449
        "INDEL1": ("1", "449", "A",     "AACGT"),
        # INDEL2: deletion [TCGA/-], anchor at ref pos 325 (unique sequence, unambiguous alignment)
        "INDEL2": ("1", "325", "GTCGA", "G"),
        # INDEL3: insertion [-/TTCC] with SNP=[D/I] (ab_a=D), anchor at ref pos 649
        "INDEL3": ("1", "649", "C",     "CTTCC"),
    }
    for name, (chrom, pos, ref, alt) in expected.items():
        assert name in by_name, f"{name} missing from position output"
        row = by_name[name]
        assert row["chromosome"] == chrom, f"{name} chromosome mismatch"
        assert row["position"] == pos, f"{name} position mismatch: got {row['position']}"
        assert row["VCF_REF"] == ref, f"{name} VCF_REF mismatch: got {row['VCF_REF']!r}"
        assert row["VCF_ALT"] == alt, f"{name} VCF_ALT mismatch: got {row['VCF_ALT']!r}"
        assert row["determination_type"], f"{name} missing determination_type"


def test_conversion_matches_expected(pipeline_output):
    rows = _read_data_rows(pipeline_output / "manifest.reference.conversion.csv")
    by_key = {(r["marker_name"], r["AB"]): r for r in rows}

    expected = {
        # SNPs
        ("SNP1", "A"): {"TOP": "A", "FORWARD": "A", "DESIGN": "A", "PLUS": "A", "VCF": "REF"},
        ("SNP1", "B"): {"TOP": "G", "FORWARD": "G", "DESIGN": "G", "PLUS": "G", "VCF": "ALT"},
        ("SNP2", "A"): {"TOP": "A", "FORWARD": "T", "DESIGN": "T", "PLUS": "T", "VCF": "REF"},
        ("SNP2", "B"): {"TOP": "C", "FORWARD": "G", "DESIGN": "G", "PLUS": "G", "VCF": "ALT"},
        ("SNP3", "A"): {"TOP": "A", "FORWARD": "T", "DESIGN": "T", "PLUS": "T", "VCF": "ALT"},
        ("SNP3", "B"): {"TOP": "C", "FORWARD": "G", "DESIGN": "G", "PLUS": "G", "VCF": "REF"},
        ("SNP4", "A"): {"TOP": "A", "FORWARD": "T", "DESIGN": "T", "PLUS": "A", "VCF": "REF"},
        ("SNP4", "B"): {"TOP": "G", "FORWARD": "C", "DESIGN": "C", "PLUS": "G", "VCF": "ALT"},
        ("SNP5", "A"): {"TOP": "A", "FORWARD": "T", "DESIGN": "A", "PLUS": "T", "VCF": "REF"},
        ("SNP5", "B"): {"TOP": "C", "FORWARD": "G", "DESIGN": "C", "PLUS": "G", "VCF": "ALT"},
        # Indels — all format fields carry I/D; VCF depends on which allele matches reference
        # INDEL1: SNP=[I/D], insertion [-/ACGT] — reference has no insertion → A(I)=ALT, B(D)=REF
        ("INDEL1", "A"): {"TOP": "I", "FORWARD": "I", "DESIGN": "I", "PLUS": "I", "VCF": "ALT"},
        ("INDEL1", "B"): {"TOP": "D", "FORWARD": "D", "DESIGN": "D", "PLUS": "D", "VCF": "REF"},
        # INDEL2: SNP=[I/D], deletion [TCGA/-] — reference has TCGA → A(I/sequence)=REF, B(D)=ALT
        ("INDEL2", "A"): {"TOP": "I", "FORWARD": "I", "DESIGN": "I", "PLUS": "I", "VCF": "REF"},
        ("INDEL2", "B"): {"TOP": "D", "FORWARD": "D", "DESIGN": "D", "PLUS": "D", "VCF": "ALT"},
        # INDEL3: SNP=[D/I], insertion [-/TTCC] — reference has no insertion → A(D)=REF, B(I)=ALT
        ("INDEL3", "A"): {"TOP": "D", "FORWARD": "D", "DESIGN": "D", "PLUS": "D", "VCF": "REF"},
        ("INDEL3", "B"): {"TOP": "I", "FORWARD": "I", "DESIGN": "I", "PLUS": "I", "VCF": "ALT"},
    }

    for key, vals in expected.items():
        assert key in by_key, f"{key} missing from conversion output"
        row = by_key[key]
        for col, val in vals.items():
            assert row[col] == val, f"{key} {col}: expected {val!r}, got {row[col]!r}"


def test_eight_variants_in_output(pipeline_output):
    rows = _read_data_rows(pipeline_output / "manifest.reference.position.csv")
    assert len(rows) == 8  # 5 SNPs + 3 indels


def test_sixteen_rows_in_conversion(pipeline_output):
    rows = _read_data_rows(pipeline_output / "manifest.reference.conversion.csv")
    assert len(rows) == 16  # 2 rows per variant × 8 variants


def test_wide_file_exists_and_has_data(pipeline_output):
    rows = _read_data_rows(pipeline_output / "manifest.reference.wide.csv")
    assert len(rows) == 8


def test_alignment_file_exists(pipeline_output):
    aln = pipeline_output / "manifest.reference.alignment.txt"
    assert aln.exists()
    text = aln.read_text()
    for name in ["SNP1", "SNP5", "INDEL1", "INDEL2", "INDEL3"]:
        assert name in text, f"{name} missing from alignment file"
