from __future__ import annotations

import gzip
import shutil
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


DEFAULT_SEQUENCE_ROLES = (
    "assembled-molecule",
    "unlocalized-scaffold",
    "unplaced-scaffold",
)


@dataclass(frozen=True)
class ReferenceDownloadStats:
    fasta_path: str
    assembly_report_path: str
    records_written: int
    size_bytes: int


def ncbi_ftp_base(accession: str, ncbi_name: str) -> str:
    prefix = accession.split(".")[0]
    if len(prefix) < 13 or "_" not in prefix:
        raise ValueError(
            "NCBI assembly accession should look like GCF_002263795.3 "
            "or GCA_000003055.4"
        )
    parts = [prefix[:3], prefix[4:7], prefix[7:10], prefix[10:13]]
    accession_dir = f"{accession}_{ncbi_name}"
    return "https://ftp.ncbi.nlm.nih.gov/genomes/all/" + "/".join(
        parts + [accession_dir]
    )


def selected_sequence_accessions(
    report_path: Path,
    sequence_roles: Iterable[str] = DEFAULT_SEQUENCE_ROLES,
) -> set[str]:
    roles = set(sequence_roles)
    accessions: set[str] = set()
    with report_path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            if fields[1] not in roles:
                continue
            accessions.add(fields[4])
            accessions.add(fields[6])
    accessions.discard("na")
    if not accessions:
        raise ValueError(f"No selected sequence accessions found in {report_path}")
    return accessions


def fasta_accession(header: str) -> str:
    first = header[1:].split(maxsplit=1)[0]
    return first.split("|")[-1]


def filter_fasta_by_accession(input_gz: Path, output_path: Path, keep: set[str]) -> int:
    written = 0
    keep_record = False
    output_path.parent.mkdir(parents=True, exist_ok=True)
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
        raise ValueError(f"No FASTA records were retained for {output_path}")
    return written


def download_ncbi_reference(
    *,
    species: str,
    assembly: str,
    accession: str,
    ncbi_name: str,
    source_root: str,
    sequence_roles: Iterable[str] = DEFAULT_SEQUENCE_ROLES,
    force: bool = False,
    downloader: Callable[[str, Path], None] | None = None,
) -> ReferenceDownloadStats:
    dest_dir = Path(source_root) / species / "references" / assembly
    fasta_path = dest_dir / f"{assembly}.refseq.fa"
    report_path = dest_dir / f"{assembly}.assembly_report.txt"
    if not force and (fasta_path.exists() or report_path.exists()):
        raise FileExistsError(
            f"Reference output already exists in {dest_dir}. Use force=True to replace it."
        )

    dest_dir.mkdir(parents=True, exist_ok=True)
    base = ncbi_ftp_base(accession, ncbi_name)
    basename = f"{accession}_{ncbi_name}"
    download = downloader or _download_url

    with tempfile.TemporaryDirectory(prefix="genotype_converter_refs_") as temp_name:
        temp_dir = Path(temp_name)
        temp_report = temp_dir / f"{basename}_assembly_report.txt"
        fasta_gz = temp_dir / f"{basename}_genomic.fna.gz"
        download(f"{base}/{basename}_assembly_report.txt", temp_report)
        download(f"{base}/{basename}_genomic.fna.gz", fasta_gz)
        shutil.copyfile(temp_report, report_path)
        keep = selected_sequence_accessions(report_path, sequence_roles)
        records_written = filter_fasta_by_accession(fasta_gz, fasta_path, keep)

    return ReferenceDownloadStats(
        fasta_path=str(fasta_path),
        assembly_report_path=str(report_path),
        records_written=records_written,
        size_bytes=fasta_path.stat().st_size,
    )


def _download_url(url: str, path: Path) -> None:
    with urllib.request.urlopen(url) as response, path.open("wb") as handle:
        shutil.copyfileobj(response, handle)
