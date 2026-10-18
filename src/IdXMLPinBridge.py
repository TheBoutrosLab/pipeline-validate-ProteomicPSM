#!/usr/bin/env python3
"""Bridge OpenMS idXML PSMs and the MSBooster/Percolator PIN workflow."""

import argparse
import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pyopenms as oms


BASE_FEATURES = [
    "SpecId", "Label", "ScanNr", "retentiontime", "deltamass",
    "ExpMass", "CalcMass", "mass", "peplen", "enzN", "enzC",
    "enzInt", "dm", "absdm", "score",
]


def load_idxml(path):
    proteins, peptides = [], oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(path), proteins, peptides)
    return proteins, peptides


def scan_number(reference):
    matches = re.findall(r"(?:^|\s)scan=(\d+)(?=\s|$)", reference)
    if len(matches) != 1:
        raise ValueError(f"Spectrum reference needs exactly one scan=N: {reference!r}")
    return int(matches[0])


def mzml_scans(path):
    scans = set()
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.rsplit("}", 1)[-1] == "spectrum":
            identifier = element.get("id", "")
            if any(child.get("accession") == "MS:1000511" and child.get("value") == "2"
                   for child in element.iter() if child.tag.rsplit("}", 1)[-1] == "cvParam"):
                scans.add(scan_number(identifier))
            element.clear()
    return scans


def as_text(value):
    return value.decode() if isinstance(value, bytes) else str(value)


def meta_text(pid, key):
    return as_text(pid.getMetaValue(key)) if pid.metaValueExists(key) else ""


def match_mzml(pid, scan, paths, scan_sets, source_stems):
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
        if index < len(source_stems):
            matches = [candidate for candidate in candidates
                       if paths[candidate].stem.lower() == source_stems[index]]
            if len(matches) == 1:
                return matches[0]
        if index in candidates:
            return index
    raise ValueError(f"Scan {scan} occurs in multiple mzML files without an unambiguous file origin")


def export(args):
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
                raise ValueError(f"PSM {pid_index}:{hit_index} lacks protein or target/decoy annotation")
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
            extra.update(part.strip() for part in as_text(search.getMetaValue("extra_features")).split(",") if part.strip())
    features = BASE_FEATURES[:9] + [f"charge{z}" for z in range(min(charges), max(charges) + 1)] + BASE_FEATURES[9:] + sorted(extra - set(BASE_FEATURES)) + ["Peptide", "Proteins"]
    oms.PercolatorInfile.store("openms.pin", peptides, [feature.encode() for feature in features],
                              args.enzyme.encode(), min(charges), max(charges))
    with open("openms.pin", newline="") as source:
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
        raise ValueError(f"OpenMS exported {len(rows)} PIN rows for {len(locations)} PSMs; cannot preserve the idXML mapping")
    include_file_name = "FileName" in fields
    fields = [field for field in fields if field != "FileName"]
    fields = fields[:-1] + ["log10_evalue", fields[-1]]
    mapping = []
    rows_by_run = {index: [] for index in range(len(paths))}
    for rank, (row, (pid_index, hit_index, scan, charge, score, higher, run_index, file_identifier)) in enumerate(zip(rows, locations), 1):
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
        with open(f"run_{run_index + 1:02d}.pin", "w", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(run_rows)
    with open(args.mapping, "w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(("SpecId", "pid_index", "hit_index"))
        writer.writerows(mapping)


def integrate(args):
    proteins, peptides = load_idxml(args.idxml)
    with open(args.mapping, newline="") as source:
        mapping = {row["SpecId"]: (int(row["pid_index"]), int(row["hit_index"]))
                   for row in csv.DictReader(source, delimiter="\t")}
    results = {}
    for filename in (args.target, args.decoy):
        with open(filename, newline="") as source:
            for row in csv.DictReader(source, delimiter="\t"):
                identifier = row["PSMId"]
                if identifier not in mapping:
                    raise ValueError(f"Percolator returned unknown PSMId {identifier!r}")
                if identifier in results:
                    raise ValueError(f"Percolator returned duplicate PSMId {identifier!r}")
                results[identifier] = (float(row["score"]), float(row["q-value"]), float(row["posterior_error_prob"]))
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
    for pid_index in range(len(peptides)):
        if pid_index not in by_pid:
            pid = peptides[pid_index]
            pid.setScoreType("Posterior Error Probability")
            pid.setHigherScoreBetter(False)
            peptides[pid_index] = pid
    oms.IdXMLFile().store(str(args.output), proteins, peptides)


def main():
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
