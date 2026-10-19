#!/usr/bin/env python3
"""Annotate idXML identifications with ProtLoom PSM and protein-group scores."""

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path

import pyopenms as oms


def string(value):
    """Decode pyOpenMS byte strings while preserving other values as text."""
    return value.decode() if isinstance(value, bytes) else str(value)


def rows(path, required):
    """Read a TSV and verify its required columns."""
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        return list(reader)


def key(native_id, sequence, charge):
    """Build the key used to match ProtLoom and idXML PSMs."""
    return native_id, sequence, int(charge)


def protloom_score(row, mode):
    """Select the ProtLoom PEP or q-value for a primary score."""
    return 1.0 - float(row["probability"]) if mode == "pep" else float(row["q_value"])


def run(args):
    """Annotate idXML peptide and protein identifications from ProtLoom tables."""
    groups = {}
    for row in rows(args.groups, ("group_id", "leading_protein", "protein_members",
                                  "protein_tier", "q_value", "probability")):
        members = [member.strip() for member in row["protein_members"].split(";") if member.strip()]
        if not members or row["leading_protein"] not in members:
            raise ValueError(f"Invalid members for ProtLoom group {row['group_id']}")
        groups[row["group_id"]] = (row, members)

    assignments = {}
    for row in rows(args.assignments, ("unimod_sequence", "charge", "group_id",
                                       "leading_protein")):
        species = (row["unimod_sequence"], int(row["charge"]))
        if species in assignments:
            raise ValueError(f"Duplicate ProtLoom assignment for {species}")
        assignments[species] = row

    psms = defaultdict(list)
    for row in rows(args.psms, ("native_id", "unimod_sequence", "charge", "rt", "mz",
                                 "protein_tier", "peptide_class", "q_value", "probability")):
        species = (row["unimod_sequence"], int(row["charge"]))
        assignment = assignments.get(species)
        if assignment is None:
            continue
        group_id = assignment["group_id"]
        if group_id not in groups:
            continue
        group, _ = groups[group_id]
        if assignment["leading_protein"] != group["leading_protein"]:
            raise ValueError(f"Leading protein disagrees for group {group_id}")
        psms[key(row["native_id"], *species)].append((row, group_id))

    proteins = []
    peptides = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(args.idxml), proteins, peptides)
    if not proteins:
        raise ValueError("Input idXML has no ProteinIdentification")
    runs = {string(protein.getIdentifier()): protein for protein in proteins}
    if len(runs) != len(proteins):
        raise ValueError("Input idXML has duplicate ProteinIdentification identifiers")
    protloom_run_ids = {}
    identifiers = set(runs)
    for run_id in runs:
        candidate = f"{run_id}_ProtLoom"
        suffix = 2
        while candidate in identifiers:
            candidate = f"{run_id}_ProtLoom_{suffix}"
            suffix += 1
        protloom_run_ids[run_id] = candidate
        identifiers.add(candidate)
    retained = oms.PeptideIdentificationList()
    matched_only = oms.PeptideIdentificationList()
    used_groups = defaultdict(set)
    unmatched_runs = set()
    matched_hits = 0
    unmatched_hits = 0

    for peptide in peptides:
        original = oms.PeptideIdentification(peptide) if args.unmatched_psms == "keep" else None
        original_hits = original.getHits() if original is not None else []
        native_id = (string(peptide.getMetaValue("spectrum_reference"))
                     if peptide.metaValueExists("spectrum_reference") else "")
        run_id = string(peptide.getIdentifier())
        if run_id not in runs:
            raise ValueError(f"PeptideIdentification references unknown run {run_id}")
        original_score_type = string(peptide.getScoreType())
        selected = []
        unmatched = []
        for hit_index, hit in enumerate(peptide.getHits()):
            sequence = string(hit.getSequence().toUniModString())
            candidates = psms.get(key(native_id, sequence, hit.getCharge()), [])
            matches = [index for index, (row, _) in enumerate(candidates)
                       if abs(float(row["rt"]) - peptide.getRT()) < 0.001
                       and abs(float(row["mz"]) - peptide.getMZ()) < 0.0001]
            if len(matches) > 1:
                raise ValueError(f"Ambiguous ProtLoom PSM for {native_id}, {sequence}")
            if not matches:
                if original is not None:
                    unmatched.append(original_hits[hit_index])
                unmatched_hits += 1
                continue
            row, group_id = candidates.pop(matches[0])
            group, _ = groups[group_id]
            leading = group["leading_protein"]
            hit.setMetaValue("protein_tier", row["protein_tier"])
            hit.setMetaValue("peptide_class", row["peptide_class"])
            hit.setMetaValue("q_value", float(row["q_value"]))
            hit.setMetaValue("probability", float(row["probability"]))
            hit.setMetaValue("protloom_pep", 1.0 - float(row["probability"]))
            hit.setMetaValue("protloom_group_id", group_id)
            if args.primary_score != "original":
                hit.setMetaValue(original_score_type or "original_primary_score", hit.getScore())
                hit.setScore(protloom_score(row, args.primary_score))
            # Existing coordinates may belong to a different protein. Keep them only
            # when they already describe the selected leading accession.
            evidence = next((item for item in hit.getPeptideEvidences()
                             if string(item.getProteinAccession()) == leading), None)
            if evidence is None:
                evidence = oms.PeptideEvidence()
                evidence.setProteinAccession(leading)
            hit.setPeptideEvidences([evidence])
            selected.append(hit)
            used_groups[run_id].add(group_id)
            matched_hits += 1
        if selected:
            annotated = (oms.PeptideIdentification(peptide)
                         if args.unmatched_psms == "keep" else peptide)
            annotated.setHits(selected)
            if args.unmatched_psms == "keep":
                annotated.setIdentifier(protloom_run_ids[run_id])
            if args.primary_score != "original":
                annotated.setScoreType(
                    "ProtLoom PEP" if args.primary_score == "pep" else "ProtLoom q-value"
                )
                annotated.setHigherScoreBetter(False)
                annotated.setSignificanceThreshold(1.0)
            retained.push_back(annotated)
            matched_only.push_back(annotated)
        if unmatched and args.unmatched_psms == "keep":
            original.setHits(unmatched)
            retained.push_back(original)
            unmatched_runs.add(run_id)

    unmatched = sum(map(len, psms.values()))
    if unmatched:
        raise ValueError(f"{unmatched} ProtLoom PSM rows could not be matched to idXML hits")

    updated_proteins = []
    matched_proteins = []
    for run_id, protein in runs.items():
        if args.unmatched_psms == "keep" and run_id in unmatched_runs:
            updated_proteins.append(protein)
        if args.unmatched_psms == "keep" and run_id not in used_groups:
            continue
        old_hits = {string(hit.getAccession()): hit for hit in protein.getHits()}
        updated = oms.ProteinIdentification()
        updated.setIdentifier(protloom_run_ids[run_id] if args.unmatched_psms == "keep" else run_id)
        updated.setDateTime(protein.getDateTime())
        updated.setSearchEngine(protein.getSearchEngine())
        updated.setSearchEngineVersion(protein.getSearchEngineVersion())
        updated.setSearchParameters(protein.getSearchParameters())
        paths = []
        protein.getPrimaryMSRunPath(paths)
        if paths:
            updated.setPrimaryMSRunPath(paths)
        meta_keys = []
        protein.getKeys(meta_keys)
        for meta_key in meta_keys:
            updated.setMetaValue(meta_key, protein.getMetaValue(meta_key))
        original_protein_score_type = string(protein.getScoreType())
        new_hits = []
        for group_id in sorted(used_groups[run_id]):
            group, members = groups[group_id]
            for accession in members:
                hit = (oms.ProteinHit(old_hits[accession])
                       if accession in old_hits else oms.ProteinHit())
                if args.primary_score != "original" and accession in old_hits:
                    hit.setMetaValue(
                        original_protein_score_type or "original_primary_score", hit.getScore()
                    )
                hit.setAccession(accession)
                hit.setScore(protloom_score(group, args.primary_score))
                hit.setMetaValue("q_value", float(group["q_value"]))
                hit.setMetaValue("probability", float(group["probability"]))
                hit.setMetaValue("protloom_pep", 1.0 - float(group["probability"]))
                hit.setMetaValue("protein_tier", group["protein_tier"])
                hit.setMetaValue("protloom_group_id", group_id)
                hit.setMetaValue("protloom_leading_protein", group["leading_protein"])
                new_hits.append(hit)
            protein_group = oms.ProteinGroup()
            # pyOpenMS 3.5 binds ProteinGroup.accessions as vector<bytes>.
            protein_group.accessions = [accession.encode("utf-8") for accession in members]
            protein_group.probability = float(group["probability"])
            updated.insertProteinGroup(protein_group)
        updated.setHits(new_hits)
        updated.setScoreType("ProtLoom protein group PEP" if args.primary_score == "pep"
                             else "ProtLoom protein group q-value")
        updated.setHigherScoreBetter(False)
        updated.setSignificanceThreshold(1.0)
        updated_proteins.append(updated)
        matched_proteins.append(updated)

    oms.IdXMLFile().store(str(args.output), updated_proteins, retained)
    if args.unmatched_psms == "keep":
        oms.IdXMLFile().store(str(args.matched_output), matched_proteins, matched_only)
    else:
        shutil.copyfile(args.output, args.matched_output)
    group_count = sum(map(len, used_groups.values()))
    print(f"Annotated {matched_hits} PSMs and {group_count} protein groups; "
          f"{unmatched_hits} unmatched PSMs "
          f"{'kept' if args.unmatched_psms == 'keep' else 'removed'}")


def main():
    """Parse command-line options and annotate the requested idXML file."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("idxml", "psms", "groups", "assignments", "output"):
        parser.add_argument(f"--{name}", required=True, type=Path)
    parser.add_argument("--matched-output", required=True, type=Path)
    parser.add_argument("--primary-score", choices=("original", "pep", "q_value", "q-value"),
                        default="original")
    parser.add_argument("--unmatched-psms", choices=("remove", "keep"), default="remove")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
