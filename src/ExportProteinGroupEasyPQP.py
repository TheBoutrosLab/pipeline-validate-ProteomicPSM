#!/usr/bin/env python3
"""Export EasyPQP tables by joining q-value PSMs to inferred protein groups.

The per-search-engine q-value idXML provides PSM scores and spectrum coordinates.
The inference output provides accepted protein groups. For Epifany, that same
output also provides the peptide-to-protein links, including any links changed
by greedy group resolution; the q-value idXML does not contain those changes.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping, MutableMapping, Sequence

import pyopenms as oms

from ProteomicIO import as_text, load_idxml, match_raw_file, meta_text, scan_number, write_rows


GroupKey = tuple[str, ...]


@dataclass(frozen=True)
class EpifanyEvidenceKey:
    """Identify a PSM within its mzML run, independent of idXML merge indices."""

    raw_name: str
    spectrum_reference: str
    modified_sequence: str
    charge: int

    @classmethod
    def from_psm(
        cls, pid: oms.PeptideIdentification, hit: oms.PeptideHit, raw_name: str
    ) -> EpifanyEvidenceKey:
        """Build the same lookup key from either Epifany or q-value idXML."""
        return cls(raw_name, meta_text(pid, "spectrum_reference"),
                   as_text(hit.getSequence().toUniModString()), hit.getCharge())


@dataclass(frozen=True)
class EpifanyEvidence:
    """Protein links and coordinates of one Epifany peptide identification."""

    rt: float
    mz: float
    accessions: frozenset[str]


@dataclass
class EpifanyEvidenceIndex:
    """Find Epifany's protein links for PSMs in the per-engine q-value idXMLs."""

    by_psm: dict[EpifanyEvidenceKey, list[EpifanyEvidence]] = field(default_factory=dict)

    @classmethod
    def from_peptides(
        cls,
        peptides: Iterable[oms.PeptideIdentification],
        raw_paths: Sequence[Path],
        source_stems: Sequence[str],
        scan_cache: MutableMapping[Path, set[int]],
    ) -> EpifanyEvidenceIndex:
        """Index Epifany's PSM-level protein links by resolved mzML origin.

        Epifany rewrites id_merge_index when merging input idXML files, so
        resolve it through Epifany's own primary MS run paths before matching.
        """
        index = cls()
        for pid in peptides:
            hits = pid.getHits()
            if not hits:
                continue
            hit = hits[0]
            reference = meta_text(pid, "spectrum_reference")
            raw = match_raw_file(pid, scan_number(reference), raw_paths,
                                 source_stems, scan_cache)
            key = EpifanyEvidenceKey.from_psm(pid, hit, raw.name)
            evidence = EpifanyEvidence(
                pid.getRT(), pid.getMZ(),
                frozenset(as_text(item.getProteinAccession())
                          for item in hit.getPeptideEvidences()),
            )
            index.by_psm.setdefault(key, []).append(evidence)
        return index

    def match(
        self, pid: oms.PeptideIdentification, hit: oms.PeptideHit, raw_name: str
    ) -> frozenset[str] | None:
        """Return matching protein links, or None if Epifany lacks this PSM.

        RT and m/z distinguish records sharing the same lookup key. Conflicting
        matches fail rather than assigning a PSM to an arbitrary group.
        """
        key = EpifanyEvidenceKey.from_psm(pid, hit, raw_name)
        matches = {item.accessions for item in self.by_psm.get(key, ())
                   if abs(item.rt - pid.getRT()) < 0.001
                   and abs(item.mz - pid.getMZ()) < 0.0001}
        if len(matches) > 1:
            raise ValueError(f"Conflicting Epifany peptide evidence for {key}")
        return next(iter(matches)) if matches else None


@dataclass
class AcceptedProteinGroups:
    """Groups passing the protein cutoff, indexed by member accession.

    ``accession_groups`` maps each protein to every accepted group containing it.
    A group key is its sorted member tuple, while ``leaders`` retains the first
    accession reported by the inference tool for EasyPQP's Protein ID field.
    ``epifany_evidence`` is present for Epifany in every resolution mode.
    """

    accession_groups: dict[str, set[GroupKey]]
    leaders: dict[GroupKey, str]
    count: int
    epifany_evidence: EpifanyEvidenceIndex | None

    @classmethod
    def load(
        cls, path: str | Path, cutoff: float, decoy_string: str, algorithm: str,
        raw_paths: Sequence[Path], scan_cache: MutableMapping[Path, set[int]],
    ) -> AcceptedProteinGroups:
        """Load target groups above the probability cutoff from inference output."""
        epifany_evidence = None
        if algorithm == "Fido":
            raw_groups = json.loads(Path(path).read_text(encoding="utf-8"))
            groups = [(float(group["probability"]),
                       [str(accession) for accession in group["accessions"]])
                      for group in raw_groups]
        else:
            runs, peptides = load_idxml(path)
            if not runs:
                raise ValueError(f"{algorithm} idXML contains no protein identifications")
            if algorithm == "Epifany":
                # Greedy resolution changes these peptide evidences in Epifany's
                # output idXML; the separate q-value idXML retains the input links.
                source_paths = []
                runs[0].getPrimaryMSRunPath(source_paths)
                source_stems = [Path(as_text(source)).stem.lower() for source in source_paths]
                epifany_evidence = EpifanyEvidenceIndex.from_peptides(
                    peptides, raw_paths, source_stems, scan_cache
                )
            groups = []
            for run in runs:
                # ProteinProphet writes protein groups; Epifany writes its
                # inferred groups as indistinguishable-protein entries.
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
            # Normalize group identity without changing the reported leader.
            identity = tuple(sorted(members))
            leaders[identity] = members[0]
            for accession in members:
                accession_groups[accession].add(identity)
            accepted += 1
        if not accepted:
            raise ValueError(f"No target {algorithm} groups pass the probability cutoff")
        return cls(accession_groups, leaders, accepted, epifany_evidence)

    def match_psm_epifany_evidence(
        self, pid: oms.PeptideIdentification, hit: oms.PeptideHit, raw_name: str
    ) -> frozenset[str] | None:
        """Return the proteins Epifany associates with this q-value PSM."""
        if self.epifany_evidence is None:
            raise ValueError("Epifany peptide evidence was not loaded")
        return self.epifany_evidence.match(pid, hit, raw_name)

    def groups_for_accessions(self, accessions: Iterable[str]) -> set[GroupKey]:
        """Find every accepted group linked to any supplied protein accession."""
        return set().union(*(self.accession_groups.get(accession, set())
                             for accession in accessions))


def export(
    manifest: Iterable[Mapping[str, Any]],
    group_path: str | Path,
    threshold: float,
    group_cutoff: float,
    decoy_string: str,
    prefix: str,
    algorithm: str,
) -> None:
    """Write cohort EasyPQP tables and counts for each exclusion stage.

    Use per-engine q-value idXML for PSM filtering and spectrum mapping. Use
    Epifany's output idXML for peptide-to-protein links when Epifany is selected;
    the other algorithms use links already present in the q-value idXML.
    """
    samples = list(manifest)
    raw_owner = {}
    all_raw_paths = []
    for sample in samples:
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
            if raw not in all_raw_paths:
                all_raw_paths.append(raw)
    selected = {}
    tied_spectra = set()
    passing_psms = 0
    unmapped_psms = 0
    scan_cache = {}
    mapped_psms = 0
    group_match_counts = [0, 0, 0]
    unmatched_epifany_psms = 0
    accepted_groups = AcceptedProteinGroups.load(
        group_path, group_cutoff, decoy_string, algorithm, all_raw_paths, scan_cache
    )

    for sample in samples:
        raw_paths = [Path(name) for name in sample["mzmls"]]

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
            native_id = meta_text(pid, "spectrum_reference")
            scan = scan_number(native_id)
            raw = match_raw_file(pid, scan, raw_paths, source_stems, scan_cache)
            # Only Epifany can revise peptide-to-protein evidence during
            # inference. Its output replaces links from the q-value idXML.
            evidence_accessions = {
                as_text(evidence.getProteinAccession())
                for evidence in hit.getPeptideEvidences()
            }
            if algorithm == "Epifany":
                resolved = accepted_groups.match_psm_epifany_evidence(pid, hit, raw.name)
                if resolved is None:
                    unmatched_epifany_psms += 1
                    evidence_accessions = set()
                else:
                    evidence_accessions = set(resolved)
            groups = accepted_groups.groups_for_accessions(evidence_accessions)
            # Count group matches before collapsing spectra or excluding
            # ambiguous peptide sequences, so resolution's effect is visible.
            group_match_counts[min(len(groups), 2)] += 1
            if len(groups) != 1:
                unmapped_psms += 1
                continue
            mapped_psms += 1
            group = next(iter(groups))
            charge = hit.getCharge()
            peptide = as_text(hit.getSequence().toUnmodifiedString())
            key = (raw.name, scan, charge)
            mapped_accessions = frozenset(
                accession for accession in evidence_accessions
                if accession and not accession.startswith(decoy_string)
            )
            candidate = (q_value, peptide, group, accepted_groups.leaders[group], raw.stem,
                         mapped_accessions)
            # A spectrum contributes one PSM: keep its best q-value, but drop
            # equal-score ties that disagree on peptide or protein group.
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
    # EasyPQP receives one group per unmodified peptide sequence. Exclude a
    # sequence if its surviving spectra point to different accepted groups.
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
            leader = accepted_groups.leaders[group]
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
        ("accepted_protein_groups", accepted_groups.count),
        ("psms_at_or_below_qvalue", passing_psms),
        ("psms_matching_zero_groups", group_match_counts[0]),
        ("psms_matching_one_group", group_match_counts[1]),
        ("psms_matching_multiple_groups", group_match_counts[2]),
        # This is the unique-group count before spectrum and peptide collapsing.
        ("psms_retained_after_resolution", mapped_psms),
        ("psms_without_epifany_match", unmatched_epifany_psms),
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


def main() -> None:
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
