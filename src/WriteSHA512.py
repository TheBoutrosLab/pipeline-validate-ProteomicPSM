#!/usr/bin/env python3
"""Write a streaming SHA-512 checksum beside an idXML file."""

import argparse
import hashlib
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("idxml", type=Path)
    args = parser.parse_args()

    digest = hashlib.sha512()
    with args.idxml.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)

    output = args.idxml.with_name(args.idxml.name + ".sha512")
    output.write_text(f"{digest.hexdigest()}  {args.idxml.name}\n", encoding="ascii")


if __name__ == "__main__":
    main()
