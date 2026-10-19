#!/usr/bin/env python3
"""Export the minimal FragPipe-style reports consumed by EasyPQP library mode."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

from IdXMLPinBridge import as_text, load_idxml, meta_text, mzml_scans, scan_number


def read_rows(path, required):
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def write_rows(path, fields, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def psm_key(native_id, unimod_sequence, charge):
    return native_id, unimod_sequence, int(charge)


def match_raw_file(pid, scan, paths, source_stems, scan_cache):
    if len(paths) == 1:
        return paths[0]
    for key in ("base_name", "file_origin"):
        origin = meta_text(pid, key)
        if origin:
            matches = [path for path in paths if path.stem.lower() == Path(origin).stem.lower()]
            if len(matches) == 1:
                return matches[0]
    merge_index = meta_text(pid, "id_merge_index")
    if merge_index:
        index = int(merge_index)
        if index < len(source_stems):
            matches = [path for path in paths if path.stem.lower() == source_stems[index]]
            if len(matches) == 1:
                return matches[0]
    candidates = []
    for path in paths:
        if path not in scan_cache:
            scan_cache[path] = mzml_scans(path)
        if scan in scan_cache[path]:
            candidates.append(path)
    if len(candidates) != 1:
        raise ValueError(
            f"Scan {scan} cannot be assigned to one mzML file from {paths}; "
            "check the idXML file-origin metadata and mzML CSV"
        )
    return candidates[0]


def export(manifest, threshold, output_prefix):
    assignments = defaultdict(set)
    selected = []
    accepted_psms = 0
    seen_raw_names = set()

    for sample in manifest:
        sample_name = sample["sample_name"]
        paths = [Path(name) for name in sample["mzmls"]]
        if not paths:
            raise ValueError(f"Sample {sample_name} has no mzML files")
        for path in paths:
            if "." in path.stem:
                raise ValueError(
                    f"EasyPQP splits Spectrum on dots; mzML run name cannot contain dots: {path.name}"
                )
            if path.name in seen_raw_names:
                raise ValueError(f"mzML basename is repeated across samples: {path.name}")
            seen_raw_names.add(path.name)

        groups = {}
        for row in read_rows(sample["group_fdr"], ("group_id", "leading_protein", "q_value")):
            if float(row["q_value"]) <= threshold:
                groups[row["group_id"]] = row["leading_protein"]

        by_species = {}
        for row in read_rows(sample["assignments"],
                             ("unmodified_sequence", "unimod_sequence", "charge",
                              "group_id", "leading_protein", "is_contam")):
            key = (row["unimod_sequence"], int(row["charge"]))
            if key in by_species:
                raise ValueError(f"Duplicate ProtLoom assignment for {sample_name}: {key}")
            by_species[key] = row

        fdr_rows = defaultdict(list)
        for row in read_rows(sample["psm_fdr"],
                             ("native_id", "unmodified_sequence", "unimod_sequence",
                              "charge", "rt", "mz", "q_value")):
            if float(row["q_value"]) > threshold:
                continue
            species = (row["unimod_sequence"], int(row["charge"]))
            assignment = by_species.get(species)
            if assignment is None:
                raise ValueError(f"No ProtLoom assignment for {sample_name}: {species}")
            if assignment["unmodified_sequence"] != row["unmodified_sequence"]:
                raise ValueError(f"ProtLoom peptide sequences disagree for {sample_name}: {species}")
            group_id = assignment["group_id"]
            if group_id not in groups:
                continue
            leading = assignment["leading_protein"]
            if leading != groups[group_id]:
                raise ValueError(f"ProtLoom leading protein differs for group {group_id}")
            key = psm_key(row["native_id"], row["unimod_sequence"], row["charge"])
            fdr_rows[key].append((row, group_id, leading))

        proteins, peptides = load_idxml(sample["idxml"])
        source_stems = []
        if proteins:
            source_paths = []
            proteins[0].getPrimaryMSRunPath(source_paths)
            source_stems = [Path(as_text(path)).stem.lower() for path in source_paths]
        scan_cache = {}
        for pid in peptides:
            hits = pid.getHits()
            if not hits:
                continue
            hit = hits[0]  # FragPipe-SpecLib's idXML reader uses the first hit.
            native_id = meta_text(pid, "spectrum_reference")
            sequence = as_text(hit.getSequence().toUniModString())
            key = psm_key(native_id, sequence, hit.getCharge())
            candidates = fdr_rows.get(key)
            if not candidates:
                continue
            if len(candidates) > 1:
                matching = [index for index, (candidate, _, _) in enumerate(candidates)
                            if abs(float(candidate["rt"]) - pid.getRT()) < 0.001
                            and abs(float(candidate["mz"]) - pid.getMZ()) < 0.0001]
                if len(matching) != 1:
                    raise ValueError(f"Cannot distinguish duplicate PSMs in {sample_name}: {key}")
                candidate_index = matching[0]
            else:
                candidate_index = 0
            row, group_id, leading = candidates.pop(candidate_index)
            scan = scan_number(native_id)
            raw_file = match_raw_file(pid, scan, paths, source_stems, scan_cache)
            peptide = as_text(hit.getSequence().toUnmodifiedString())
            if peptide != row["unmodified_sequence"]:
                raise ValueError(f"idXML and ProtLoom PSM differ for {sample_name}: {key}")
            assignments[peptide].add((group_id, leading))
            selected.append((peptide, group_id, leading, raw_file, scan, hit.getCharge()))
            accepted_psms += 1

        unmatched = sum(len(values) for values in fdr_rows.values())
        if unmatched:
            raise ValueError(
                f"{unmatched} accepted ProtLoom PSMs in {sample_name} were not found as "
                "first hits in its idXML; EasyPQP cannot read them from that file"
            )

    ambiguous = {peptide for peptide, groups in assignments.items() if len(groups) > 1}
    psm_rows = []
    peptide_rows = []
    for peptide, group_id, leading, raw_file, scan, charge in selected:
        if peptide in ambiguous:
            continue
        psm_rows.append({
            "Spectrum": f"{raw_file.stem}.{scan}.{scan}.{charge}",
            "Spectrum File": raw_file.name,
            "Peptide": peptide,
        })
    for peptide in sorted(assignments):
        if peptide not in ambiguous:
            _, leading = next(iter(assignments[peptide]))
            peptide_rows.append({
                "Peptide": peptide,
                # ProtLoom does not expose a gene in these tables. A stable
                # accession keeps EasyPQP's required Gene field nonempty.
                "Gene": leading,
                "Protein ID": leading,
            })

    if not psm_rows:
        raise ValueError("No PSMs remain after ProtLoom FDR and ambiguity filtering")

    psm_rows.sort(key=lambda row: (row["Spectrum File"], row["Spectrum"], row["Peptide"]))
    write_rows(f"{output_prefix}-psm.tsv", ("Spectrum", "Spectrum File", "Peptide"), psm_rows)
    write_rows(f"{output_prefix}-peptides.tsv", ("Peptide", "Gene", "Protein ID"), peptide_rows)
    summary = (
        ("accepted_psms_before_ambiguity", accepted_psms),
        ("ambiguous_peptide_sequences", len(ambiguous)),
        ("psms_removed_for_ambiguity", accepted_psms - len(psm_rows)),
        ("exported_psms", len(psm_rows)),
        ("exported_peptide_sequences", len(peptide_rows)),
    )
    write_rows(f"{output_prefix}-easypqp_export_summary.tsv", ("metric", "count"),
               ({"metric": name, "count": value} for name, value in summary))
    for name, value in summary:
        print(f"{name}: {value}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--fdr", type=float, default=0.01)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()
    if not 0 <= args.fdr <= 1:
        parser.error("--fdr must be between 0 and 1")
    if not args.output_prefix or args.output_prefix in (".", "..") or any(
        character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
        for character in args.output_prefix
    ):
        parser.error("--output-prefix must be a valid dataset ID")
    export(json.loads(args.manifest.read_text()), args.fdr, args.output_prefix)


if __name__ == "__main__":
    main()
