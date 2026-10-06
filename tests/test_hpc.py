from diflow.hpc import write_slurm_script


def test_write_slurm_script(tmp_path):
    path = write_slurm_script(
        tmp_path / "run.slurm",
        command="diflow infer --help",
        job_name="DIFLOW test",
        cpus=4,
        mem_gb=16,
        hours=12,
        email="user@example.com",
    )
    text = path.read_text(encoding="utf-8")
    assert "#SBATCH --job-name=DIFLOW_test" in text
    assert "#SBATCH --cpus-per-task=4" in text
    assert "#SBATCH --mem=16G" in text
    assert "#SBATCH --time=12:00:00" in text
    assert "diflow infer --help" in text
    assert "set -euo pipefail" in text
