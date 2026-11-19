#!/usr/bin/env python3
"""Run one NFTest case from the validation pipeline repository root."""

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
import sys

import yaml
from dotenv import load_dotenv


CLUSTERS = {"", "boutroslab_aws", "pines"}


def resolve_run_root(repo: Path, cluster: str) -> Path:
    """Resolve the test base directory and add the AWS Slurm job directory."""
    root = Path(os.environ.get("NFT_RUN_ROOT") or str(repo)).expanduser()
    root = (repo / root).resolve()
    if cluster == "boutroslab_aws":
        job_id = os.environ.get("SLURM_JOB_ID")
        if not job_id:
            raise ValueError("boutroslab_aws requires SLURM_JOB_ID; run inside a Slurm job")
        root /= job_id
    return root


def prepare_config(repo: Path, test_case: str, test_dir: Path, cluster: str) -> Path:
    """Validate the case and select Apptainer for Pines."""
    config_path = repo / "nftest.yml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    matches = [case for case in config["cases"] if case["name"] == test_case]
    if len(matches) != 1:
        raise ValueError(f"Expected one NFTest case named {test_case!r}; found {len(matches)}")
    test_dir.mkdir(parents=True, exist_ok=True)
    if cluster == "pines":
        profiles = list(matches[0].get("profiles") or [])
        if "apptainer" not in profiles:
            profiles.append("apptainer")
        matches[0]["profiles"] = profiles
        config_path = test_dir / "nftest.yml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    return config_path


def run_test(repo: Path, test_case: str, root: Path, cluster: str) -> int:
    """Run NFTest from the active Python environment."""
    case_id = hashlib.sha256(test_case.encode("utf-8")).hexdigest()[:12]
    test_dir = root / ".nftest-temp" / case_id
    config_path = prepare_config(repo, test_case, test_dir, cluster)
    env = {
        **os.environ,
        "NFT_CLUSTER": cluster,
        "NFT_INIT": str(repo),
        "NFT_PIPELINE": str(repo),
        "NFT_TESTDIR": str(repo),
        "NFT_OUTPUT": str(root),
        "NFT_TEMP": str(test_dir),
    }
    print(f"Repository:  {repo}", flush=True)
    print(f"Test case:   {test_case}", flush=True)
    print(f"Cluster:     {cluster or 'default'}", flush=True)
    print(f"Work dir:    {test_dir}", flush=True)
    print(f"Output root: {root}", flush=True)
    result = subprocess.run(
        [sys.executable, "-m", "nftest", "run", "-c", str(config_path), test_case],
        cwd=test_dir,
        env=env,
        check=False,
    )
    return result.returncode if result.returncode >= 0 else 128 - result.returncode


def main() -> int:
    """Load test settings and report configuration errors."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("test_case")
    args = parser.parse_args()
    repo = Path.cwd()
    if not (repo / "nftest.yml").is_file():
        parser.error("launch from the validation pipeline repository root")
    load_dotenv(repo / ".env", override=True)
    cluster = os.environ.get("NFT_CLUSTER", "").strip()
    try:
        if cluster not in CLUSTERS:
            raise ValueError(f"Unknown NFT_CLUSTER {cluster!r}; use boutroslab_aws or pines")
        root = resolve_run_root(repo, cluster)
        return run_test(repo, args.test_case, root, cluster)
    except (ValueError, OSError, yaml.YAMLError) as error:
        parser.print_usage(file=sys.stderr)
        print(f"{parser.prog}: error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
