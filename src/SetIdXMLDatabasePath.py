#!/usr/bin/env python3
"""Use one database path in every idXML before cohort-wide Prophet analysis."""

import argparse
from pathlib import Path

import pyopenms as oms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    proteins = []
    peptides = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(str(args.input), proteins, peptides)
    if not proteins:
        raise ValueError(f"{args.input} has no ProteinIdentification")
    database = str(args.database.resolve(strict=True))
    for protein in proteins:
        settings = protein.getSearchParameters()
        settings.db = database
        protein.setSearchParameters(settings)
    oms.IdXMLFile().store(str(args.output), proteins, peptides)


if __name__ == "__main__":
    main()
