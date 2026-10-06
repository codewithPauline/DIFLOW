import gzip
from pathlib import Path

import pandas as pd
import pytest

from diflow.io import allele_counts_from_vcf, read_popmap


VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\tA2\tB1\tB2
1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1
1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t./.\t0/0\t0/1
1\t30\t.\tG\tA,C\t.\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/0
"""


def test_read_popmap_and_count_alleles(tmp_path: Path):
    popmap_path = tmp_path / "popmap.tsv"
    popmap_path.write_text(
        "sample\tpopulation\nA1\tA\nA2\tA\nB1\tB\nB2\tB\n",
        encoding="utf-8",
    )
    vcf_path = tmp_path / "tiny.vcf"
    vcf_path.write_text(VCF, encoding="utf-8")

    popmap = read_popmap(popmap_path)
    counts = allele_counts_from_vcf(vcf_path, popmap)

    locus10 = counts[counts["pos"] == 10].set_index("population")
    assert locus10.loc["A", "alt_count"] == 1
    assert locus10.loc["A", "called_chromosomes"] == 4
    assert locus10.loc["B", "alt_count"] == 3
    assert locus10.loc["B", "called_chromosomes"] == 4

    locus20 = counts[counts["pos"] == 20].set_index("population")
    assert locus20.loc["A", "called_chromosomes"] == 2
    assert locus20.loc["B", "called_chromosomes"] == 4

    assert 30 not in counts["pos"].tolist()


def test_popmap_duplicate_rejected(tmp_path: Path):
    path = tmp_path / "bad.tsv"
    path.write_text("sample\tpopulation\nX\tA\nX\tB\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_popmap(path)



def test_gzip_vcf_input(tmp_path):
    vcf = tmp_path / "tiny.vcf.gz"
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n"
        "1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1/1\n"
    )
    with gzip.open(vcf, "wt", encoding="utf-8") as handle:
        handle.write(text)

    popmap = pd.DataFrame(
        {"sample": ["s1", "s2"], "population": ["A", "B"]}
    )
    counts = allele_counts_from_vcf(vcf, popmap)
    assert len(counts) == 2
    assert counts["called_chromosomes"].tolist() == [2, 2]
