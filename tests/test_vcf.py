from __future__ import annotations

from genotype_converter.vcf import write_site_vcf


def test_build_writes_site_only_vcf_with_nucleotide_alleles(pipeline_output):
    vcf_path = pipeline_output / "manifest.reference.sites.vcf"
    assert vcf_path.is_file()

    lines = vcf_path.read_text().splitlines()
    header = next(line for line in lines if line.startswith("#CHROM"))
    assert header.split("\t") == [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
    ]

    records = [line.split("\t") for line in lines if not line.startswith("#")]
    assert records
    assert all(len(record) == 8 for record in records)
    snp1 = next(record for record in records if record[2] == "SNP1")
    assert snp1[:5] == ["1", "300", "SNP1", "A", "G"]


def test_site_vcf_normalizes_multi_alt_nucleotide_alleles(tmp_path):
    output = tmp_path / "sites.vcf"
    stats = write_site_vcf(
        [
            {
                "marker_name": "SNP_MULTI",
                "alt_marker_name": "",
                "chromosome": "2",
                "position": "10",
                "ref_allele": "C",
                "alt_allele": "A/G",
                "determination_type": "SNP_ALIGNED",
                "A_in_PLUS": "A",
                "B_in_PLUS": "G",
                "A_vcf": "ALT",
                "B_vcf": "ALT",
            },
            {
                "marker_name": "UNPOSITIONED",
                "chromosome": "",
                "position": "",
                "ref_allele": "",
                "alt_allele": "",
            },
        ],
        str(output),
    )

    assert stats.records_total == 2
    assert stats.records_written == 1
    assert stats.records_skipped == 1
    record = next(line for line in output.read_text().splitlines() if not line.startswith("#"))
    assert record.split("\t")[:5] == ["2", "10", "SNP_MULTI", "C", "A,G"]
    assert "A_ALLELE=A" in record
    assert "B_ALLELE=G" in record
