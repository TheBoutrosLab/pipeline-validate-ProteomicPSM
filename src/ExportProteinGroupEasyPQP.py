#!/usr/bin/env python3
"""Export EasyPQP tables from OpenMS PSM q-values and inferred protein groups."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

from ExportEasyPQPReports import match_raw_file, write_rows
from IdXMLPinBridge import as_text, load_idxml, meta_text, scan_number


def protein_groups(path, cutoff, decoy_string, algorithm):
    """Index accepted target protein groups by their member accessions."""
    if algorithm == "Fido":
        raw_groups = json.loads(Path(path).read_text())
        groups = [(float(group["probability"]),
                   [str(accession) for accession in group["accessions"]])
                  for group in raw_groups]
    else:
        runs, _ = load_idxml(path)
        if not runs:
            raise ValueError(f"{algorithm} idXML contains no protein identifications")
        groups = []
        for run in runs:
            run_groups = (run.getProteinGroups() if algorithm == "ProteinProphet"
                          else run.getIndistinguishableProteins())
            groups.extend((group.probability,
                           [as_text(accession) for accession in group.accessions])
                          for group in run_groups)
    accession_groups = defaultdict(set)
    leaders = {}
    accepted = 0
    for probability, members in groups:
        if probability < cutoff:
            continue
        if not members or any(member.startswith(decoy_string) for member in members):
            continue
        identity = tuple(sorted(members))
        leaders[identity] = members[0]
        for accession in members:
            accession_groups[accession].add(identity)
        accepted += 1
    if not accepted:
        raise ValueError(f"No target {algorithm} groups pass the probability cutoff")
    return accession_groups, leaders, accepted


def export(manifest, group_path, threshold, group_cutoff, decoy_string, prefix,
           algorithm):
    """Write PSM, peptide, and exclusion-summary tables for a cohort."""
    accession_groups, leaders, accepted_groups = protein_groups(
        group_path, group_cutoff, decoy_string, algorithm
    )
    selected = {}
    tied_spectra = set()
    passing_psms = 0
    unmapped_psms = 0
    raw_owner = {}
    scan_cache = {}
    mapped_psms = 0

    for sample in manifest:
        sample_name = sample["sample_name"]
        raw_paths = [Path(name) for name in sample["mzmls"]]
        if not raw_paths:
            raise ValueError(f"Sample {sample_name} has no mzML files")
        for raw in raw_paths:
            if "." in raw.stem:
                raise ValueError(f"mzML run name cannot contain dots for EasyPQP: {raw.name}")
            owner = raw_owner.setdefault(raw.name, sample_name)
            if owner != sample_name:
                raise ValueError(f"mzML basename occurs in multiple samples: {raw.name}")

        runs, peptides = load_idxml(sample["idxml"])
        source_stems = []
        if runs:
            source_paths = []
            runs[0].getPrimaryMSRunPath(source_paths)
            source_stems = [Path(as_text(path)).stem.lower() for path in source_paths]
        for pid in peptides:
            hits = pid.getHits()
            if not hits:
                continue
            hit = hits[0]
            score_type = as_text(pid.getScoreType()).lower()
            if pid.isHigherScoreBetter() or not any(
                marker in score_type for marker in ("q-value", "q_value", "q value")
            ):
                raise ValueError(f"Expected lower-is-better q-values in {sample['idxml']}")
            q_value = float(hit.getScore())
            if q_value > threshold:
                continue
            if hit.metaValueExists("target_decoy") and "decoy" in as_text(
                hit.getMetaValue("target_decoy")
            ).lower():
                continue
            passing_psms += 1
            evidence_accessions = {
                as_text(evidence.getProteinAccession())
                for evidence in hit.getPeptideEvidences()
            }
            groups = set().union(*(accession_groups.get(accession, set())
                                   for accession in evidence_accessions))
            if len(groups) != 1:
                unmapped_psms += 1
                continue
            mapped_psms += 1
            group = next(iter(groups))
            native_id = meta_text(pid, "spectrum_reference")
            scan = scan_number(native_id)
            raw = match_raw_file(pid, scan, raw_paths, source_stems, scan_cache)
            charge = hit.getCharge()
            peptide = as_text(hit.getSequence().toUnmodifiedString())
            key = (raw.name, scan, charge)
            mapped_accessions = frozenset(
                accession for accession in evidence_accessions
                if accession and not accession.startswith(decoy_string)
            )
            candidate = (q_value, peptide, group, leaders[group], raw.stem,
                         mapped_accessions)
            previous = selected.get(key)
            if previous is None or q_value < previous[0]:
                selected[key] = candidate
                tied_spectra.discard(key)
            elif q_value == previous[0] and (peptide, group) != (previous[1], previous[2]):
                tied_spectra.add(key)
            elif q_value == previous[0]:
                selected[key] = (*previous[:5], previous[5] | mapped_accessions)

    for key in tied_spectra:
        selected.pop(key, None)
    peptide_groups = defaultdict(set)
    peptide_accessions = defaultdict(set)
    for _, peptide, group, _, _, mapped_accessions in selected.values():
        peptide_groups[peptide].add(group)
        peptide_accessions[peptide].update(mapped_accessions)
    ambiguous = {peptide for peptide, groups in peptide_groups.items() if len(groups) > 1}

    psm_rows = []
    peptide_rows = []
    for (raw_name, scan, charge), (_, peptide, group, _, raw_stem, _) in sorted(selected.items()):
        if peptide in ambiguous:
            continue
        psm_rows.append({"Spectrum": f"{raw_stem}.{scan}.{scan}.{charge}",
                         "Spectrum File": raw_name, "Peptide": peptide})
    for peptide in sorted(peptide_groups):
        if peptide not in ambiguous:
            group = next(iter(peptide_groups[peptide]))
            leader = leaders[group]
            peptide_rows.append({"Peptide": peptide, "Gene": leader,
                                 "Protein ID": leader,
                                 "Mapped Proteins": ";".join(sorted(
                                     peptide_accessions[peptide] - {leader}
                                 ))})
    if not psm_rows:
        raise ValueError("No PSMs remain after q-value and protein-group filtering")

    write_rows(f"{prefix}-psm.tsv", ("Spectrum", "Spectrum File", "Peptide"), psm_rows)
    write_rows(f"{prefix}-peptides.tsv",
               ("Peptide", "Gene", "Protein ID", "Mapped Proteins"), peptide_rows)
    summary = [
        ("accepted_protein_groups", accepted_groups),
        ("psms_at_or_below_qvalue", passing_psms),
        ("psms_without_unique_group", unmapped_psms),
        ("conflicting_tied_spectra", len(tied_spectra)),
        ("duplicate_spectrum_rows_collapsed", mapped_psms - len(tied_spectra) -
                                              len(selected)),
        ("ambiguous_peptide_sequences", len(ambiguous)),
        ("psms_removed_for_ambiguity", len(selected) - len(psm_rows)),
        ("exported_psms", len(psm_rows)),
        ("exported_peptide_sequences", len(peptide_rows)),
    ]
    write_rows(f"{prefix}-easypqp_export_summary.tsv", ("metric", "count"),
               ({"metric": metric, "count": count} for metric, count in summary))
    for metric, count in summary:
        print(f"{metric}: {count}")


def main():
    """Parse command-line options and export the requested cohort tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--groups", required=True, type=Path)
    parser.add_argument("--algorithm", required=True,
                        choices=("ProteinProphet", "Fido", "Epifany"))
    parser.add_argument("--psm-qvalue", required=True, type=float)
    parser.add_argument("--min-protein-probability", required=True, type=float)
    parser.add_argument("--decoy-string", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()
    if not 0 <= args.psm_qvalue <= 1:
        parser.error("--psm-qvalue must be between zero and one")
    if not 0 <= args.min_protein_probability < 1:
        parser.error("--min-protein-probability must be in [0, 1)")
    if not args.output_prefix or any(c not in
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
        for c in args.output_prefix):
        parser.error("--output-prefix must be a valid dataset ID")
    export(json.loads(args.manifest.read_text()), args.groups, args.psm_qvalue,
           args.min_protein_probability, args.decoy_string, args.output_prefix,
           args.algorithm)


if __name__ == "__main__":
    main()
