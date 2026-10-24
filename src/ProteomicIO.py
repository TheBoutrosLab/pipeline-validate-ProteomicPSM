#!/usr/bin/env python3
"""Shared idXML, mzML, and TSV helpers for proteomic report exporters."""

from __future__ import annotations

import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Iterable, Mapping, MutableMapping, Sequence

import pyopenms as oms


def as_text(value: object) -> str:
    """Decode pyOpenMS byte strings while preserving other values as text."""
    return value.decode() if isinstance(value, bytes) else str(value)


def load_idxml(
    path: str | Path,
) -> tuple[list[oms.ProteinIdentification], oms.PeptideIdentificationList]:
    """Load protein and peptide identifications from an idXML file."""
    proteins, peptides = [], oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(path), proteins, peptides)
    return proteins, peptides


def meta_text(pid: oms.PeptideIdentification, key: str) -> str:
    """Read an optional peptide identification metadata value as text."""
    return as_text(pid.getMetaValue(key)) if pid.metaValueExists(key) else ""


def scan_number(reference: str) -> int:
    """Extract exactly one scan number from a spectrum reference."""
    matches = re.findall(r"(?:^|\s)scan=(\d+)(?=\s|$)", reference)
    if len(matches) != 1:
        raise ValueError(f"Spectrum reference needs exactly one scan=N: {reference!r}")
    return int(matches[0])


def mzml_scans(path: str | Path) -> set[int]:
    """Collect MS2 scan numbers from an mzML file."""
    scans = set()
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.rsplit("}", 1)[-1] == "spectrum":
            identifier = element.get("id", "")
            if any(child.get("accession") == "MS:1000511" and child.get("value") == "2"
                   for child in element.iter() if child.tag.rsplit("}", 1)[-1] == "cvParam"):
                scans.add(scan_number(identifier))
            element.clear()
    return scans


def write_rows(
    path: str | Path, fields: Sequence[str], rows: Iterable[Mapping[str, object]]
) -> None:
    """Write rows to a UTF-8 TSV file."""
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def match_raw_file(
    pid: oms.PeptideIdentification,
    scan: int,
    paths: Sequence[Path],
    source_stems: Sequence[str],
    scan_cache: MutableMapping[Path, set[int]],
) -> Path:
    """Find a PSM's mzML file from provenance or an unambiguous scan number."""
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
