#!/usr/bin/env python3
"""Bridge OpenMS idXML PSMs and the MSBooster/Percolator PIN workflow."""

import argparse
import csv
import math
from pathlib import Path

import pyopenms as oms

from ProteomicIO import as_text, load_idxml, meta_text, mzml_ms2_native_ids, mzml_scans, scan_number


BASE_FEATURES = [
    "SpecId", "Label", "ScanNr", "retentiontime", "deltamass",
    "ExpMass", "CalcMass", "mass", "peplen", "enzN", "enzC",
    "enzInt", "dm", "absdm", "score",
]


def msbooster_score(engine, score_type, higher_is_better, score):
    """Return the PIN score column and value expected by MSBooster for an engine."""
    engine = engine.lower()
    if not math.isfinite(score):
        raise ValueError(f"{engine} has a non-finite PSM score: {score}")
    if engine in ("msfragger", "comet"):
        if score_type != "expect" or higher_is_better:
            raise ValueError(
                f"{engine} requires lower-is-better 'expect' idXML scores; "
                f"got {score_type!r} (higher_is_better={higher_is_better})"
            )
        if score < 0:
            raise ValueError(f"{engine} has a negative expectation value: {score}")
        # An expectation value of zero is rounded underflow; keep its log finite.
        return "log10_evalue", math.log10(max(score, 1e-30))
    if engine == "sage":
        if score_type != "ln(hyperscore)" or not higher_is_better:
            raise ValueError(
                "Sage requires higher-is-better 'ln(hyperscore)' idXML scores; "
                f"got {score_type!r} (higher_is_better={higher_is_better})"
            )
        try:
            hyperscore = math.exp(score)
        except OverflowError as error:
            raise ValueError(f"Sage hyperscore overflows: {score}") from error
        if not math.isfinite(hyperscore) or hyperscore <= 0:
            raise ValueError(f"Sage has an invalid hyperscore: {score}")
        return "hyperscore", hyperscore
    raise ValueError(f"MSBooster PIN export does not support search engine {engine!r}")


def match_mzml(pid, scan, paths, scan_sets, native_id_sets, source_stems):
    """Find the raw file by native ID when available, otherwise by scan number."""
    native_id = meta_text(pid, "pepxml_spectrum_name")
    if native_id:
        candidates = [index for index, ids in enumerate(native_id_sets) if native_id in ids]
        if not candidates:
            raise ValueError(
                f"Native spectrum ID {native_id!r} from pepXML is absent from the sample's mzML files"
            )
    else:
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
    identifier = f"native spectrum ID {native_id!r}" if native_id else f"scan {scan}"
    raise ValueError(
        f"{identifier} occurs in multiple mzML files without an unambiguous file origin"
    )


def export(args):
    """Export idXML hits to per-run PIN files and record their original positions."""
    proteins, peptides = load_idxml(args.idxml)
    if not peptides:
        raise ValueError("Input idXML has no peptide identifications")
    paths = [Path(path) for path in args.mzml]
    if len({path.name.lower() for path in paths}) != len(paths):
        raise ValueError("mzML files for one sample must have distinct basenames")
    native_id_sets = [mzml_ms2_native_ids(path) for path in paths]
    scan_sets = [mzml_scans(path) for path in paths] if any(
        not meta_text(pid, "pepxml_spectrum_name") for pid in peptides
    ) else [set() for _ in paths]
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
        run_index = match_mzml(pid, scan, paths, scan_sets, native_id_sets, source_stems)
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
            score_column, score = msbooster_score(
                args.search_engine, as_text(pid.getScoreType()),
                pid.isHigherScoreBetter(), hit.getScore()
            )
            locations.append((pid_index, hit_index, scan, charge, score_column,
                              score, run_index, file_identifier))
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
    expected_column = "hyperscore" if args.search_engine.lower() == "sage" else "log10_evalue"
    if expected_column not in fields:
        fields = fields[:-1] + [expected_column, fields[-1]]
    mapping = []
    rows_by_run = {index: [] for index in range(len(paths))}
    for rank, (row, location) in enumerate(zip(rows, locations), 1):
        (pid_index, hit_index, scan, charge, score_column, score,
         run_index, file_identifier) = location
        if int(row["ScanNr"]) != scan:
            raise ValueError("OpenMS PIN order or scan numbers differ from the idXML")
        if include_file_name and row.pop("FileName") != file_identifier:
            raise ValueError("OpenMS PIN source-file order differs from the idXML")
        spec_id = f"run_{run_index + 1:02d}.{scan}.{scan}.{charge}_{rank}"
        row["SpecId"] = spec_id
        row[score_column] = str(score)
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
    exp.add_argument("--search-engine", required=True)
    exp.add_argument("--mzml", nargs="+", required=True)
    imp = sub.add_parser("integrate")
    for name in ("idxml", "mapping", "target", "decoy", "output"):
        imp.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    {"export": export, "integrate": integrate}[args.command](args)


if __name__ == "__main__":
    main()
