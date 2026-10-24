#!/usr/bin/env python3
"""Bridge OpenMS idXML PSMs and the MSBooster/Percolator PIN workflow."""

import argparse
import csv
from pathlib import Path

import pyopenms as oms

from ProteomicIO import as_text, load_idxml, meta_text, mzml_scans, scan_number


BASE_FEATURES = [
    "SpecId", "Label", "ScanNr", "retentiontime", "deltamass",
    "ExpMass", "CalcMass", "mass", "peplen", "enzN", "enzC",
    "enzInt", "dm", "absdm", "score",
]


def match_mzml(pid, scan, paths, scan_sets, source_stems):
    """Find the raw file for a scan, requiring provenance if scans repeat."""
    candidates = [index for index, scans in enumerate(scan_sets) if scan in scans]
    if not candidates:
        raise ValueError(f"Scan {scan} in idXML is absent from the sample's mzML files")
    if len(candidates) == 1:
        return candidates[0]
    origin = meta_text(pid, "file_origin")
    if origin:
        matches = [index for index in candidates
                   if paths[index].stem.lower() == Path(origin).stem.lower()]
        if len(matches) == 1:
            return matches[0]
    merge_index = meta_text(pid, "id_merge_index")
    if merge_index:
        index = int(merge_index)
        if 0 <= index < len(source_stems):
            matches = [candidate for candidate in candidates
                       if paths[candidate].stem.lower() == source_stems[index]]
            if len(matches) == 1:
                return matches[0]
    raise ValueError(
        f"Scan {scan} occurs in multiple mzML files without an unambiguous file origin"
    )


def export(args):
    """Export idXML hits to per-run PIN files and record their original positions."""
    proteins, peptides = load_idxml(args.idxml)
    if not peptides:
        raise ValueError("Input idXML has no peptide identifications")
    paths = [Path(path) for path in args.mzml]
    if len({path.name.lower() for path in paths}) != len(paths):
        raise ValueError("mzML files for one sample must have distinct basenames")
    scan_sets = [mzml_scans(path) for path in paths]
    source_stems = []
    if proteins:
        source_paths = []
        proteins[0].getPrimaryMSRunPath(source_paths)
        source_stems = [Path(as_text(path)).stem.lower() for path in source_paths]
    for run_index, path in enumerate(paths, 1):
        alias = Path(f"run_{run_index:02d}.mzML")
        if alias.exists() and alias.resolve() != path.resolve():
            raise ValueError(f"Cannot create mzML alias {alias}: name already exists")
        if not alias.exists():
            alias.symlink_to(path.resolve())
    locations = []
    charges = []
    for pid_index, pid in enumerate(peptides):
        scan = scan_number(meta_text(pid, "spectrum_reference"))
        run_index = match_mzml(pid, scan, paths, scan_sets, source_stems)
        file_identifier = meta_text(pid, "file_origin")
        merge_index = meta_text(pid, "id_merge_index")
        if merge_index:
            file_identifier += f"|{merge_index}"
        for hit_index, hit in enumerate(pid.getHits()):
            if not hit.getPeptideEvidences() or not hit.metaValueExists("target_decoy"):
                raise ValueError(
                    f"PSM {pid_index}:{hit_index} lacks protein or target/decoy annotation"
                )
            charge = hit.getCharge()
            if charge < 1:
                raise ValueError(f"PSM {pid_index}:{hit_index} has invalid charge {charge}")
            locations.append((pid_index, hit_index, scan, charge, hit.getScore(),
                              pid.isHigherScoreBetter(), run_index, file_identifier))
            charges.append(charge)
    if not locations:
        raise ValueError("Input idXML contains no PSMs")
    extra = set()
    for protein in proteins:
        search = protein.getSearchParameters()
        if search.metaValueExists("extra_features"):
            extra.update(
                part.strip()
                for part in as_text(search.getMetaValue("extra_features")).split(",")
                if part.strip()
            )
    features = (BASE_FEATURES[:9]
                + [f"charge{z}" for z in range(min(charges), max(charges) + 1)]
                + BASE_FEATURES[9:] + sorted(extra - set(BASE_FEATURES))
                + ["Peptide", "Proteins"])
    oms.PercolatorInfile.store("openms.pin", peptides, [feature.encode() for feature in features],
                              args.enzyme.encode(), min(charges), max(charges))
    with open("openms.pin", newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source, delimiter="\t")
        rows = list(reader)
        fields = reader.fieldnames
    if not fields or fields[-1] != "Proteins":
        raise ValueError("OpenMS PIN has an unexpected header")
    for row in rows:
        # OpenMS emits multiple protein accessions as extra trailing TSV cells.
        if None in row:
            row["Proteins"] = ";".join([row["Proteins"]] + row.pop(None))
    if len(rows) != len(locations):
        raise ValueError(
            f"OpenMS exported {len(rows)} PIN rows for {len(locations)} PSMs; "
            "cannot preserve the idXML mapping"
        )
    include_file_name = "FileName" in fields
    fields = [field for field in fields if field != "FileName"]
    fields = fields[:-1] + ["log10_evalue", fields[-1]]
    mapping = []
    rows_by_run = {index: [] for index in range(len(paths))}
    for rank, (row, location) in enumerate(zip(rows, locations), 1):
        (pid_index, hit_index, scan, charge, score, higher,
         run_index, file_identifier) = location
        if int(row["ScanNr"]) != scan:
            raise ValueError("OpenMS PIN order or scan numbers differ from the idXML")
        if include_file_name and row.pop("FileName") != file_identifier:
            raise ValueError("OpenMS PIN source-file order differs from the idXML")
        spec_id = f"run_{run_index + 1:02d}.{scan}.{scan}.{charge}_{rank}"
        # MSBooster uses this value to choose calibration PSMs; monotonicity is sufficient.
        row["SpecId"] = spec_id
        row["log10_evalue"] = str(max(-30.0, min(30.0, -score if higher else score)))
        row["retentiontime"] = str(float(row["retentiontime"]) / 60.0)
        mapping.append((spec_id, pid_index, hit_index))
        rows_by_run[run_index].append(row)
    for run_index, run_rows in rows_by_run.items():
        if not run_rows:
            continue
        with open(f"run_{run_index + 1:02d}.pin", "w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(run_rows)
    with open(args.mapping, "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(("SpecId", "pid_index", "hit_index"))
        writer.writerows(mapping)


def integrate(args):
    """Write Percolator scores into their original idXML hits."""
    proteins, peptides = load_idxml(args.idxml)
    with open(args.mapping, newline="", encoding="utf-8") as source:
        mapping = {row["SpecId"]: (int(row["pid_index"]), int(row["hit_index"]))
                   for row in csv.DictReader(source, delimiter="\t")}
    results = {}
    for filename in (args.target, args.decoy):
        with open(filename, newline="", encoding="utf-8") as source:
            for row in csv.DictReader(source, delimiter="\t"):
                identifier = row["PSMId"]
                if identifier not in mapping:
                    raise ValueError(f"Percolator returned unknown PSMId {identifier!r}")
                if identifier in results:
                    raise ValueError(f"Percolator returned duplicate PSMId {identifier!r}")
                results[identifier] = (
                    float(row["score"]), float(row["q-value"]),
                    float(row["posterior_error_prob"])
                )
    by_pid = {}
    for identifier, (pid_index, hit_index) in mapping.items():
        by_pid.setdefault(pid_index, []).append((identifier, hit_index))
    for pid_index, entries in by_pid.items():
        pid = peptides[pid_index]
        hits = pid.getHits()
        old_score_type = pid.getScoreType()
        for identifier, hit_index in entries:
            hit = hits[hit_index]
            svm, qvalue, pep = results.get(identifier, (-100.0, 1.0, 1.0))
            hit.setMetaValue(old_score_type, hit.getScore())
            hit.setMetaValue("MS:1001492", svm)
            hit.setMetaValue("MS:1001491", qvalue)
            hit.setMetaValue("MS:1001493", pep)
            hit.setScore(pep)
            hits[hit_index] = hit
        pid.setHits(hits)
        pid.setScoreType("Posterior Error Probability")
        pid.setHigherScoreBetter(False)
        peptides[pid_index] = pid
    for pid_index, pid in enumerate(peptides):
        if pid_index not in by_pid:
            pid.setScoreType("Posterior Error Probability")
            pid.setHigherScoreBetter(False)
            peptides[pid_index] = pid
    oms.IdXMLFile().store(str(args.output), proteins, peptides)


def main():
    """Parse command-line options and run the requested conversion."""
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    exp = sub.add_parser("export")
    for name in ("idxml", "mapping", "enzyme"):
        exp.add_argument(f"--{name}", required=True)
    exp.add_argument("--mzml", nargs="+", required=True)
    imp = sub.add_parser("integrate")
    for name in ("idxml", "mapping", "target", "decoy", "output"):
        imp.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    {"export": export, "integrate": integrate}[args.command](args)


if __name__ == "__main__":
    main()
