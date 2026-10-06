from pathlib import Path

import pandas as pd

from diflow.pipeline import run_infer_pipeline


VCF = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA1\tA2\tB1\tB2
1\t10\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\t1/1\t0/1
1\t20\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/1
"""


def test_prepare_only_pipeline(tmp_path: Path):
    vcf = tmp_path / "tiny.vcf"
    vcf.write_text(VCF, encoding="utf-8")

    popmap = tmp_path / "popmap.tsv"
    popmap.write_text(
        "sample\tpopulation\nA1\tA\nA2\tA\nB1\tB\nB2\tB\n",
        encoding="utf-8",
    )

    coords = tmp_path / "coords.csv"
    pd.DataFrame(
        {
            "population": ["A", "B"],
            "latitude": [39.0, 39.2],
            "longitude": [-84.0, -83.8],
        }
    ).to_csv(coords, index=False)

    out = tmp_path / "results"
    result = run_infer_pipeline(
        vcf_path=vcf,
        popmap_path=popmap,
        coordinates_path=coords,
        output_dir=out,
        projection_chromosomes=4,
        k_nearest=1,
        prepare_only=True,
    )

    assert len(result.candidate_pairs) == 1
    assert result.pairwise_results.iloc[0]["status"] == "prepared"
    assert (out / "allele_counts.csv").exists()
    assert (out / "candidate_pairs.csv").exists()
    assert (out / "pairwise_results.csv").exists()
    assert (out / "run_metadata.csv").exists()
    assert len(list((out / "spectra").glob("*.npy"))) == 1
