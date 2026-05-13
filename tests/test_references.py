from __future__ import annotations

import gzip

import pytest

from genotype_converter.references import (
    filter_fasta_by_accession,
    fasta_accession,
    ncbi_ftp_base,
    selected_sequence_accessions,
)


def test_ncbi_ftp_base_uses_accession_directory_layout():
    url = ncbi_ftp_base("GCF_002263795.3", "ARS-UCD2.0")
    assert url == (
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
        "GCF/002/263/795/GCF_002263795.3_ARS-UCD2.0"
    )


def test_selected_sequence_accessions_keeps_default_roles(tmp_path):
    report = tmp_path / "assembly_report.txt"
    report.write_text(
        "\n".join(
            [
                "# comment",
                "1\tassembled-molecule\tna\tna\tCM000001.1\tna\tNC_000001.1\t=\t1\tChromosome",
                "2\talt-scaffold\tna\tna\tNW_ALT.1\tna\tNW_ALT_REF.1\t=\tna\tna",
                "3\tunplaced-scaffold\tna\tna\tNW_KEEP.1\tna\tNW_KEEP_REF.1\t=\tna\tna",
            ]
        )
        + "\n"
    )

    assert selected_sequence_accessions(report) == {
        "CM000001.1",
        "NC_000001.1",
        "NW_KEEP.1",
        "NW_KEEP_REF.1",
    }


def test_filter_fasta_by_accession_writes_selected_records(tmp_path):
    fasta_gz = tmp_path / "input.fna.gz"
    with gzip.open(fasta_gz, "wt") as handle:
        handle.write(">NC_000001.1 chromosome\nACGT\n")
        handle.write(">NW_SKIP.1 scaffold\nTTTT\n")
        handle.write(">gi|1|ref|NW_KEEP.1 scaffold\nGGGG\n")

    output = tmp_path / "filtered.fa"
    written = filter_fasta_by_accession(
        fasta_gz,
        output,
        {"NC_000001.1", "NW_KEEP.1"},
    )

    assert written == 2
    assert output.read_text() == (
        ">NC_000001.1 chromosome\nACGT\n"
        ">gi|1|ref|NW_KEEP.1 scaffold\nGGGG\n"
    )
    assert fasta_accession(">gi|1|ref|NW_KEEP.1 scaffold") == "NW_KEEP.1"


def test_filter_fasta_by_accession_rejects_empty_output(tmp_path):
    fasta_gz = tmp_path / "input.fna.gz"
    with gzip.open(fasta_gz, "wt") as handle:
        handle.write(">NC_000001.1 chromosome\nACGT\n")

    with pytest.raises(ValueError, match="No FASTA records"):
        filter_fasta_by_accession(fasta_gz, tmp_path / "filtered.fa", {"missing"})
