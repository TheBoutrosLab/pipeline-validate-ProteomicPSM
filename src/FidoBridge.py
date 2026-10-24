#!/usr/bin/env python3
"""Convert Percolator-scored idXML to Fido inputs and parse Fido groups."""

import argparse
import json
import math
import re
from pathlib import Path

from ProteomicIO import as_text, load_idxml


GROUP_LINE = re.compile(r"^\s*([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)\s*\{([^}]*)\}\s*$")


def peptide_probability(peptide_id, hit):
    """Return a posterior probability after checking the primary PEP score."""
    score_type = as_text(peptide_id.getScoreType()).lower()
    normalized = re.sub(r"[\s_-]", "", score_type)
    if normalized.endswith("score"):
        normalized = normalized[:-5]
    if normalized not in ("pep", "posteriorerrorprobability") or peptide_id.isHigherScoreBetter():
        raise ValueError(f"Fido requires a lower-is-better PEP score; found {score_type!r}")
    pep = float(hit.getScore())
    if not math.isfinite(pep) or not 0 <= pep <= 1:
        raise ValueError(f"Invalid PEP for Fido: {pep}")
    return 1.0 - pep


def prepare(idxml_paths, decoy_prefix, graph_path, proteins_path, mapping_path):
    """Write Fido's PSM graph, target/decoy lists, and accession mapping."""
    if not decoy_prefix:
        raise ValueError("A decoy accession prefix is required for Fido")
    proteins = set()
    records = []
    for path in sorted(idxml_paths):
        runs, peptides = load_idxml(path)
        for run in runs:
            proteins.update(as_text(hit.getAccession()) for hit in run.getHits())
        for peptide_id in peptides:
            hits = peptide_id.getHits()
            if not hits:
                continue
            hit = min(hits, key=lambda item: item.getScore())
            sequence = as_text(hit.getSequence().toUnmodifiedString())
            accessions = sorted({as_text(item.getProteinAccession())
                                 for item in hit.getPeptideEvidences()})
            accessions = [accession for accession in accessions if accession]
            if not sequence or not accessions:
                continue
            probability = peptide_probability(peptide_id, hit)
            proteins.update(accessions)
            records.append((sequence, accessions, probability))
    proteins.discard("")
    if not records:
        raise ValueError("No indexed PSMs with Percolator PEP were found for Fido")
    targets = sorted(protein for protein in proteins if not protein.startswith(decoy_prefix))
    decoys = sorted(protein for protein in proteins if protein.startswith(decoy_prefix))
    if not targets or not decoys:
        raise ValueError("Fido parameter estimation requires target and decoy proteins")

    accessions = {f"P{index:08d}": accession
                  for index, accession in enumerate(sorted(proteins), 1)}
    tokens = {accession: token for token, accession in accessions.items()}
    with graph_path.open("w") as graph:
        for sequence, members, probability in records:
            graph.write(f"e {sequence}\n")
            for member in members:
                graph.write(f"r {tokens[member]}\n")
            graph.write(f"p {probability:.17g}\n")
    with proteins_path.open("w") as output:
        output.write("{ " + " , ".join(tokens[name] for name in targets) + " }\n")
        output.write("{ " + " , ".join(tokens[name] for name in decoys) + " }\n")
    mapping_path.write_text(json.dumps(accessions, sort_keys=True) + "\n")
    print(f"Fido graph: {len(records)} PSMs, {len(targets)} targets, {len(decoys)} decoys")


def parse_results(results_path, mapping_path, output_path):
    """Translate Fido's native group lines back to original accessions."""
    accessions = json.loads(mapping_path.read_text())
    groups = []
    for line in results_path.read_text().splitlines():
        match = GROUP_LINE.fullmatch(line)
        if not match:
            continue
        probability = float(match.group(1))
        if not math.isfinite(probability) or not 0 <= probability <= 1:
            raise ValueError(f"Invalid Fido group probability: {probability}")
        tokens = [token.strip() for token in match.group(2).split(",")]
        if not tokens or any(not token for token in tokens):
            continue
        try:
            members = [accessions[token] for token in tokens]
        except KeyError as error:
            raise ValueError(f"Unknown protein token in Fido output: {error}") from error
        groups.append({"probability": probability, "accessions": members})
    if not groups:
        raise ValueError("Fido produced no parseable protein groups")
    output_path.write_text(json.dumps(groups, indent=2) + "\n")
    print(f"Parsed {len(groups)} Fido protein groups")


def main():
    """Parse the requested graph-export or group-import command."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--idxml", nargs="+", required=True, type=Path)
    prepare_parser.add_argument("--decoy-prefix", required=True)
    prepare_parser.add_argument("--graph", required=True, type=Path)
    prepare_parser.add_argument("--proteins", required=True, type=Path)
    prepare_parser.add_argument("--mapping", required=True, type=Path)
    parse_parser = commands.add_parser("parse")
    parse_parser.add_argument("--results", required=True, type=Path)
    parse_parser.add_argument("--mapping", required=True, type=Path)
    parse_parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare(args.idxml, args.decoy_prefix, args.graph, args.proteins,
                args.mapping)
    else:
        parse_results(args.results, args.mapping, args.output)


if __name__ == "__main__":
    main()
