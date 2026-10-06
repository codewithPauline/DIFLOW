"""Minimal VCF genotype-to-population allele counter.

The parser intentionally supports a conservative subset of VCF needed for the
first DIFLOW real-data layer: biallelic records with a GT field. Multiallelic
records and genotypes containing allele indices greater than 1 are skipped.
"""

from __future__ import annotations

from collections import defaultdict
import gzip
from pathlib import Path

import pandas as pd


def _parse_gt(sample_field: str, gt_index: int) -> tuple[int, int] | None:
    fields = sample_field.split(":")
    if gt_index >= len(fields):
        return None
    gt = fields[gt_index]
    if gt in {".", "./.", ".|."}:
        return None

    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles):
        return None
    try:
        values = [int(a) for a in alleles]
    except ValueError:
        return None
    if any(a not in (0, 1) for a in values):
        return None

    alt = sum(values)
    called = len(values)
    return alt, called


def allele_counts_from_vcf(
    vcf_path: str | Path,
    popmap: pd.DataFrame,
    *,
    min_called_chromosomes: int = 1,
) -> pd.DataFrame:
    """Aggregate biallelic GT allele counts by population for every locus.

    Parameters
    ----------
    vcf_path
        Plain-text VCF or gzip-compressed .vcf.gz path. Indexed random access is
        not required for the current streaming reader.
    popmap
        DataFrame with sample and population columns.
    min_called_chromosomes
        Minimum number of called chromosomes required to emit a
        locus-population row.

    Returns
    -------
    pandas.DataFrame
        Long-form table with chrom, pos, ref, alt, population, ref_count,
        alt_count and called_chromosomes.
    """
    if min_called_chromosomes < 1:
        raise ValueError("min_called_chromosomes must be at least 1.")
    if not {"sample", "population"}.issubset(popmap.columns):
        raise ValueError("popmap must contain sample and population columns.")

    mapping = dict(zip(popmap["sample"].astype(str), popmap["population"].astype(str)))
    path = Path(vcf_path)
    if not path.exists():
        raise FileNotFoundError(path)

    rows: list[dict] = []
    samples: list[str] | None = None

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("##"):
                continue
            line = raw.rstrip("\n")

            if line.startswith("#CHROM"):
                header = line.split("\t")
                samples = header[9:]
                missing = sorted(set(mapping) - set(samples))
                if missing:
                    raise ValueError(
                        "popmap contains samples absent from VCF: " + ", ".join(missing)
                    )
                continue

            if not line or line.startswith("#"):
                continue
            if samples is None:
                raise ValueError("VCF header line beginning #CHROM was not found.")

            fields = line.split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt = fields[:5]
            if "," in alt:
                continue

            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                continue
            gt_index = format_fields.index("GT")

            by_pop: dict[str, list[int]] = defaultdict(lambda: [0, 0])
            for sample, sample_field in zip(samples, fields[9:]):
                population = mapping.get(sample)
                if population is None:
                    continue
                parsed = _parse_gt(sample_field, gt_index)
                if parsed is None:
                    continue
                alt_count, called = parsed
                by_pop[population][0] += alt_count
                by_pop[population][1] += called

            for population, (alt_count, called) in by_pop.items():
                if called < min_called_chromosomes:
                    continue
                rows.append(
                    {
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": alt,
                        "population": population,
                        "ref_count": called - alt_count,
                        "alt_count": alt_count,
                        "called_chromosomes": called,
                    }
                )

    return pd.DataFrame(
        rows,
        columns=[
            "chrom",
            "pos",
            "ref",
            "alt",
            "population",
            "ref_count",
            "alt_count",
            "called_chromosomes",
        ],
    )
