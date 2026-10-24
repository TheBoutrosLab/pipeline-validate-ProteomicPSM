"""Exercise EasyPQP group assignment with a shared-peptide Epifany fixture."""

import importlib.util
import json
import sys
import types
import unittest
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
FIXTURE = json.loads((ROOT / "test/fixtures/epifany_shared_peptide.json").read_text())


class Sequence:
    """Expose the sequence methods used by the exporter."""

    def __init__(self, value):
        self.value = value

    def toUniModString(self):
        return self.value

    def toUnmodifiedString(self):
        return self.value


class Hit:
    """Represent one PSM hit and its current protein associations."""

    def __init__(self, row, accessions):
        self.row = row
        self.accessions = accessions

    def getSequence(self):
        return Sequence(self.row["peptide"])

    def getCharge(self):
        return 2

    def getScore(self):
        return self.row["q_value"]

    def getPeptideEvidences(self):
        return [types.SimpleNamespace(getProteinAccession=lambda name=name: name)
                for name in self.accessions]

    def metaValueExists(self, _key):
        return False


class Peptide:
    """Represent one identification with stable spectrum coordinates."""

    def __init__(self, row, accessions):
        self.row = row
        self.hit = Hit(row, accessions)

    def getHits(self):
        return [self.hit]

    def getRT(self):
        return float(self.row["scan"])

    def getMZ(self):
        return 500.0

    def getScoreType(self):
        return "q-value"

    def isHigherScoreBetter(self):
        return False


class EpifanyExportTest(unittest.TestCase):
    """Check that resolution changes group mapping but not q-value selection."""

    @classmethod
    def setUpClass(cls):
        pyopenms = types.ModuleType("pyopenms")
        pyopenms.PeptideIdentification = Peptide
        pyopenms.PeptideHit = Hit
        io = types.ModuleType("ProteomicIO")
        for name in ("as_text", "load_idxml", "match_raw_file", "meta_text",
                     "scan_number", "write_rows"):
            setattr(io, name, lambda *_args: None)
        with patch.dict(sys.modules, {"pyopenms": pyopenms, "ProteomicIO": io}):
            spec = importlib.util.spec_from_file_location(
                "ExportProteinGroupEasyPQP_test", ROOT / "src/ExportProteinGroupEasyPQP.py"
            )
            module = importlib.util.module_from_spec(spec)
            with patch.dict(sys.modules, {spec.name: module}):
                spec.loader.exec_module(module)
        cls.exporter = module

    def run_mode(self, mode):
        """Export the same q-value PSMs with one Epifany evidence mode."""
        groups = [types.SimpleNamespace(probability=0.99, accessions=members)
                  for members in FIXTURE["groups"]]
        run = types.SimpleNamespace(
            getIndistinguishableProteins=lambda: groups,
            getPrimaryMSRunPath=lambda paths: paths.append("run.mzML"),
        )
        q_peptides = [Peptide(row, row["input"]) for row in FIXTURE["psms"]]
        inferred = [Peptide(row, row[mode]) for row in FIXTURE["psms"]]
        tables = {}

        def load(path):
            return ([run], inferred) if path == "groups.idXML" else ([run], q_peptides)

        def metadata(pid, _key):
            return f"scan={pid.row['scan']}"

        with patch.multiple(
            self.exporter,
            as_text=lambda value: str(value),
            load_idxml=load,
            meta_text=metadata,
            scan_number=lambda native: int(native.split("=")[1]),
            match_raw_file=lambda *_args: Path("run.mzML"),
            write_rows=lambda path, _columns, rows: tables.setdefault(path, list(rows)),
        ):
            self.exporter.export(
                [{"sample_name": "sample", "idxml": "q.idXML", "mzmls": ["run.mzML"]}],
                "groups.idXML", 0.01, 0.01, "DECOY_", "test", "Epifany"
            )
        summary = {row["metric"]: row["count"]
                   for row in tables["test-easypqp_export_summary.tsv"]}
        return tables["test-psm.tsv"], summary

    def test_greedy_resolution_recovers_shared_psm(self):
        """Default excludes the shared PSM; resolution assigns it to one group."""
        default_psms, default_summary = self.run_mode("none")
        resolved_psms, resolved_summary = self.run_mode("resolved")
        self.assertEqual([row["Peptide"] for row in default_psms], ["UNIQUE"])
        self.assertEqual([row["Peptide"] for row in resolved_psms],
                         ["UNIQUE", "SHARED"])
        self.assertEqual(default_summary["psms_matching_multiple_groups"], 1)
        self.assertEqual(resolved_summary["psms_matching_multiple_groups"], 0)
        self.assertEqual(resolved_summary["psms_retained_after_resolution"], 2)
        self.assertEqual(resolved_summary["psms_at_or_below_qvalue"], 2)


if __name__ == "__main__":
    unittest.main()
