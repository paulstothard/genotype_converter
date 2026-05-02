from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path

from .convert_genotypes import FORMATS, _convert_allele


PLINK_FORMATS = FORMATS - {"VCF"}


@dataclass
class PlinkConvertStats:
    variants_total: int
    variants_converted: int
    variants_missing_lookup: int
    alleles_changed: int
    genotype_path: str
    variant_path: str
    sample_path: str

    @property
    def bed_path(self) -> str:
        return self.genotype_path

    @property
    def bim_path(self) -> str:
        return self.variant_path

    @property
    def fam_path(self) -> str:
        return self.sample_path


def _prefix_path(prefix: str, suffix: str) -> Path:
    return Path(f"{prefix}.{suffix}")


def _require_plink_formats(from_fmt: str, to_fmt: str) -> None:
    unknown = [fmt for fmt in (from_fmt, to_fmt) if fmt not in PLINK_FORMATS]
    if unknown:
        supported = ", ".join(sorted(PLINK_FORMATS))
        raise ValueError(
            "PLINK conversion supports allele-label formats only. "
            f"Unsupported format(s): {', '.join(unknown)}. Supported: {supported}"
        )


def _require_bfile(prefix: str) -> tuple[Path, Path, Path]:
    bed = _prefix_path(prefix, "bed")
    bim = _prefix_path(prefix, "bim")
    fam = _prefix_path(prefix, "fam")
    missing = [str(path) for path in (bed, bim, fam) if not path.exists()]
    if missing:
        raise ValueError(
            "PLINK binary conversion requires all three files: .bed, .bim, and .fam. "
            f"Missing: {', '.join(missing)}"
        )
    return bed, bim, fam


def _require_pfile(prefix: str) -> tuple[Path, Path, Path]:
    pgen = _prefix_path(prefix, "pgen")
    pvar = _prefix_path(prefix, "pvar")
    psam = _prefix_path(prefix, "psam")
    missing = [str(path) for path in (pgen, pvar, psam) if not path.exists()]
    if missing:
        raise ValueError(
            "PLINK 2 conversion requires all three files: .pgen, .pvar, and .psam. "
            f"Missing: {', '.join(missing)}"
        )
    return pgen, pvar, psam


def convert_plink_bfile(
    input_prefix: str,
    output_prefix: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    *,
    update_position: bool = False,
    require_all_markers: bool = True,
) -> PlinkConvertStats:
    """
    Convert allele labels in a PLINK binary fileset.

    The genotype matrix in .bed is unchanged. The .fam file is copied unchanged.
    The .bim allele columns are rewritten according to the lookup table.
    """
    _require_plink_formats(from_fmt, to_fmt)
    input_bed, input_bim, input_fam = _require_bfile(input_prefix)

    output_bed = _prefix_path(output_prefix, "bed")
    output_bim = _prefix_path(output_prefix, "bim")
    output_fam = _prefix_path(output_prefix, "fam")
    output_bed.parent.mkdir(parents=True, exist_ok=True)
    output_bim.parent.mkdir(parents=True, exist_ok=True)
    output_fam.parent.mkdir(parents=True, exist_ok=True)

    variants_total = 0
    variants_converted = 0
    variants_missing_lookup = 0
    alleles_changed = 0
    missing_markers: list[str] = []

    with input_bim.open() as in_handle, output_bim.open("w") as out_handle:
        for line_number, line in enumerate(in_handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            fields = stripped.split()
            if len(fields) != 6:
                raise ValueError(
                    f"{input_bim} line {line_number} has {len(fields)} fields; expected 6"
                )

            chrom, marker, cm, pos, allele1, allele2 = fields
            variants_total += 1
            rows = table.get(marker)
            if not rows:
                variants_missing_lookup += 1
                missing_markers.append(marker)
                out_handle.write("\t".join(fields) + "\n")
                continue

            allele1_out = _convert_allele(allele1, marker, from_fmt, to_fmt, table)
            allele2_out = _convert_allele(allele2, marker, from_fmt, to_fmt, table)
            if allele1_out != allele1:
                alleles_changed += 1
            if allele2_out != allele2:
                alleles_changed += 1

            lookup_row = rows[0]
            if update_position:
                chrom = lookup_row.get("chromosome", chrom) or chrom
                pos = lookup_row.get("position", pos) or pos

            out_handle.write(
                "\t".join([chrom, marker, cm, pos, allele1_out, allele2_out]) + "\n"
            )
            variants_converted += 1

    if missing_markers and require_all_markers:
        output_bim.unlink(missing_ok=True)
        preview = ", ".join(sorted(set(missing_markers))[:10])
        more = f" and {len(set(missing_markers)) - 10} more" if len(set(missing_markers)) > 10 else ""
        raise ValueError(
            "PLINK .bim contains marker(s) not found in the lookup table: "
            f"{preview}{more}"
        )

    shutil.copyfile(input_bed, output_bed)
    shutil.copyfile(input_fam, output_fam)

    return PlinkConvertStats(
        variants_total=variants_total,
        variants_converted=variants_converted,
        variants_missing_lookup=variants_missing_lookup,
        alleles_changed=alleles_changed,
        genotype_path=str(output_bed),
        variant_path=str(output_bim),
        sample_path=str(output_fam),
    )


def _pvar_header_indexes(header: list[str]) -> dict[str, int]:
    normalized = ["CHROM" if col == "#CHROM" else col for col in header]
    required = ["CHROM", "POS", "ID", "REF", "ALT"]
    missing = [col for col in required if col not in normalized]
    if missing:
        raise ValueError(f"PLINK 2 .pvar header is missing column(s): {', '.join(missing)}")
    return {col: normalized.index(col) for col in required}


def convert_plink_pfile(
    input_prefix: str,
    output_prefix: str,
    table: dict,
    from_fmt: str,
    to_fmt: str,
    *,
    update_position: bool = False,
    require_all_markers: bool = True,
) -> PlinkConvertStats:
    """
    Convert allele labels in a PLINK 2 pfile.

    The genotype matrix in .pgen is unchanged. The .psam file is copied unchanged.
    The .pvar REF/ALT columns are rewritten according to the lookup table.
    """
    _require_plink_formats(from_fmt, to_fmt)
    input_pgen, input_pvar, input_psam = _require_pfile(input_prefix)

    output_pgen = _prefix_path(output_prefix, "pgen")
    output_pvar = _prefix_path(output_prefix, "pvar")
    output_psam = _prefix_path(output_prefix, "psam")
    output_pgen.parent.mkdir(parents=True, exist_ok=True)
    output_pvar.parent.mkdir(parents=True, exist_ok=True)
    output_psam.parent.mkdir(parents=True, exist_ok=True)

    variants_total = 0
    variants_converted = 0
    variants_missing_lookup = 0
    alleles_changed = 0
    missing_markers: list[str] = []
    indexes: dict[str, int] | None = None

    with input_pvar.open() as in_handle, output_pvar.open("w") as out_handle:
        for line_number, line in enumerate(in_handle, start=1):
            stripped = line.rstrip("\n")
            if not stripped:
                out_handle.write(line)
                continue
            if stripped.startswith("##"):
                out_handle.write(line)
                continue
            if stripped.startswith("#CHROM"):
                header = stripped.split()
                indexes = _pvar_header_indexes(header)
                out_handle.write("\t".join(header) + "\n")
                continue
            if indexes is None:
                raise ValueError(
                    f"{input_pvar} line {line_number} appears before a #CHROM header. "
                    "Only headered .pvar files are supported."
                )

            fields = stripped.split()
            max_index = max(indexes.values())
            if len(fields) <= max_index:
                raise ValueError(
                    f"{input_pvar} line {line_number} has {len(fields)} fields; "
                    f"expected at least {max_index + 1}"
                )

            marker = fields[indexes["ID"]]
            alt = fields[indexes["ALT"]]
            variants_total += 1
            if "," in alt:
                raise ValueError(
                    f"{input_pvar} line {line_number} marker {marker!r} is multiallelic; "
                    "PLINK 2 conversion currently supports biallelic records only"
                )

            rows = table.get(marker)
            if not rows:
                variants_missing_lookup += 1
                missing_markers.append(marker)
                out_handle.write("\t".join(fields) + "\n")
                continue

            ref = fields[indexes["REF"]]
            ref_out = _convert_allele(ref, marker, from_fmt, to_fmt, table)
            alt_out = _convert_allele(alt, marker, from_fmt, to_fmt, table)
            if ref_out != ref:
                alleles_changed += 1
            if alt_out != alt:
                alleles_changed += 1
            fields[indexes["REF"]] = ref_out
            fields[indexes["ALT"]] = alt_out

            lookup_row = rows[0]
            if update_position:
                fields[indexes["CHROM"]] = lookup_row.get("chromosome", fields[indexes["CHROM"]]) or fields[indexes["CHROM"]]
                fields[indexes["POS"]] = lookup_row.get("position", fields[indexes["POS"]]) or fields[indexes["POS"]]

            out_handle.write("\t".join(fields) + "\n")
            variants_converted += 1

    if indexes is None:
        output_pvar.unlink(missing_ok=True)
        raise ValueError(f"{input_pvar} is missing a #CHROM header")

    if missing_markers and require_all_markers:
        output_pvar.unlink(missing_ok=True)
        preview = ", ".join(sorted(set(missing_markers))[:10])
        more = f" and {len(set(missing_markers)) - 10} more" if len(set(missing_markers)) > 10 else ""
        raise ValueError(
            "PLINK 2 .pvar contains marker(s) not found in the lookup table: "
            f"{preview}{more}"
        )

    shutil.copyfile(input_pgen, output_pgen)
    shutil.copyfile(input_psam, output_psam)

    return PlinkConvertStats(
        variants_total=variants_total,
        variants_converted=variants_converted,
        variants_missing_lookup=variants_missing_lookup,
        alleles_changed=alleles_changed,
        genotype_path=str(output_pgen),
        variant_path=str(output_pvar),
        sample_path=str(output_psam),
    )
