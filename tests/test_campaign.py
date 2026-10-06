import json

from diflow.campaign import write_validation_campaign


def test_validation_campaign_writes_manifest_and_jobs(tmp_path):
    manifest_path = write_validation_campaign(
        tmp_path / "campaign",
        replicates=5,
        starts=2,
        bootstrap_replicates=10,
        cpus=4,
        mem_gb=16,
        hours=12,
        seed=7,
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert set(manifest["jobs"]) == {
        "recovery_grid",
        "decision_evidence",
        "linked_calibration",
        "mechanistic_linkage",
    }
    assert manifest["settings"]["replicates"] == 5
    assert "release_review_command" in manifest["postprocessing"]
    assert "release-review" in manifest["postprocessing"]["release_review_command"]

    root = manifest_path.parent
    assert (root / "01_recovery_grid.slurm").exists()
    assert (root / "02_decision_evidence.slurm").exists()
    assert (root / "03_linked_calibration.slurm").exists()
    assert (root / "04_mechanistic_linkage.slurm").exists()
    assert (root / "submit_all.sh").exists()

    submit = (root / "submit_all.sh").read_text(encoding="utf-8")
    assert "sbatch" in submit
