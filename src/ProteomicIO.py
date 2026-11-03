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
    """Extract one numeric scan or diaTracer frame ID from a spectrum reference."""
    matches = re.findall(r"(?:^|\s)(?:scan|frame)=(\d+)(?=\s|$)", reference)
    if len(matches) != 1:
        raise ValueError(f"Spectrum reference needs exactly one scan=N or frame=N: {reference!r}")
    return int(matches[0])


def mzml_ms2_native_ids(path: str | Path) -> set[str]:
    """Collect exact native IDs of MS2 spectra from an mzML file."""
    identifiers = set()
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.rsplit("}", 1)[-1] == "spectrum":
            if any(child.get("accession") == "MS:1000511" and child.get("value") == "2"
                   for child in element.iter() if child.tag.rsplit("}", 1)[-1] == "cvParam"):
                native_id = element.get("id", "")
                if not native_id or native_id in identifiers:
                    raise ValueError(f"Missing or duplicate MS2 native ID {native_id!r} in {path}")
                identifiers.add(native_id)
            element.clear()
    return identifiers


def mzml_scans(path: str | Path) -> set[int]:
    """Collect numeric MS2 scan IDs for inputs without a preserved native ID."""
    scans = set()
    for _, element in ET.iterparse(path, events=("end",)):
        if element.tag.rsplit("}", 1)[-1] == "spectrum":
            identifier = element.get("id", "")
            if any(child.get("accession") == "MS:1000511" and child.get("value") == "2"
                   for child in element.iter() if child.tag.rsplit("}", 1)[-1] == "cvParam"):
                # A frame number is not interchangeable with a scan number.
                if re.search(r"(?:^|\s)scan=\d+(?=\s|$)", identifier):
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
        if 0 <= index < len(source_stems):
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
