#!/usr/bin/env python3
from __future__ import annotations

import gzip
import shutil
import sys
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parent / "sources"


@dataclass(frozen=True)
class ReferenceDownload:
    species: str
    assembly: str
    accession: str
    ncbi_name: str

    @property
    def ftp_base(self) -> str:
        prefix = self.accession.split(".")[0]
        parts = [prefix[:3], prefix[4:7], prefix[7:10], prefix[10:13]]
        accession_dir = f"{self.accession}_{self.ncbi_name}"
        return "https://ftp.ncbi.nlm.nih.gov/genomes/all/" + "/".join(parts + [accession_dir])

    @property
    def basename(self) -> str:
        return f"{self.accession}_{self.ncbi_name}"


REFERENCES = [
    ReferenceDownload("bos_taurus", "ARS_UCD1_2", "GCF_002263795.1", "ARS-UCD1.2"),
    ReferenceDownload("bos_taurus", "ARS_UCD_v2_0", "GCF_002263795.3", "ARS-UCD2.0"),
    ReferenceDownload("bos_taurus", "UMD3_1", "GCF_000003055.4", "Bos_taurus_UMD_3.1"),
    ReferenceDownload("sus_scrofa", "Sscrofa11_1", "GCF_000003025.6", "Sscrofa11.1"),
]


def download(url: str, path: Path) -> None:
    print(f"Downloading {url}", flush=True)
    with urllib.request.urlopen(url) as response, path.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def assembled_molecule_accessions(report_path: Path) -> set[str]:
    accessions: set[str] = set()
    with report_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sequence_role = fields[1]
            genbank_accession = fields[4]
            refseq_accession = fields[6]
            if sequence_role == "assembled-molecule":
                accessions.add(genbank_accession)
                accessions.add(refseq_accession)
    accessions.discard("na")
    if not accessions:
        raise RuntimeError(f"No assembled-molecule accessions found in {report_path}")
    return accessions


def fasta_accession(header: str) -> str:
    first = header[1:].split(maxsplit=1)[0]
    return first.split("|")[-1]


def filter_fasta_by_accession(input_gz: Path, output_path: Path, keep: set[str]) -> int:
    written = 0
    keep_record = False
    with gzip.open(input_gz, "rt") as source, output_path.open("w") as out:
        for line in source:
            if line.startswith(">"):
                keep_record = fasta_accession(line) in keep
                if keep_record:
                    written += 1
                    out.write(line)
            elif keep_record:
                out.write(line)
    if written == 0:
        raise RuntimeError(f"No FASTA records were retained for {output_path}")
    return written


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="genotype_converter_refs_") as temp_name:
        temp_dir = Path(temp_name)
        for reference in REFERENCES:
            dest_dir = SOURCE_ROOT / reference.species / "references" / reference.assembly
            dest_dir.mkdir(parents=True, exist_ok=True)
            output_path = dest_dir / f"{reference.assembly}.chromosomes.fa"
            report_path = dest_dir / f"{reference.assembly}.assembly_report.txt"
            fasta_gz = temp_dir / f"{reference.basename}_genomic.fna.gz"
            temp_report = temp_dir / f"{reference.basename}_assembly_report.txt"

            download(f"{reference.ftp_base}/{reference.basename}_assembly_report.txt", temp_report)
            download(f"{reference.ftp_base}/{reference.basename}_genomic.fna.gz", fasta_gz)
            shutil.copyfile(temp_report, report_path)
            keep = assembled_molecule_accessions(report_path)
            written = filter_fasta_by_accession(fasta_gz, output_path, keep)
            size_gb = output_path.stat().st_size / (1024**3)
            print(
                f"Wrote {output_path} with {written} assembled molecule(s), "
                f"{size_gb:.2f} GiB",
                flush=True,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
