"""Extract the OpenMS digestion enzyme rule from an idXML file."""

import argparse
import json

from pyopenms import IdXMLFile, PeptideIdentificationList


def parse_args():
    """Parse the idXML input and JSON output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-path', required=True)
    parser.add_argument('--output-path', required=True)
    return parser.parse_args()


def extract_rule(input_path):
    """Return the one digestion enzyme shared by all protein identifications."""
    proteins = []
    peptides = PeptideIdentificationList()
    IdXMLFile().load(input_path, proteins, peptides)
    if not proteins:
        raise ValueError('idXML contains no ProteinIdentification.')

    enzymes = set()
    for identification in proteins:
        enzyme = identification.getSearchParameters().digestion_enzyme
        enzymes.add((enzyme.getName(), enzyme.getRegEx()))
    if len(enzymes) != 1:
        raise ValueError(f'Expected one digestion enzyme in idXML, found {enzymes!r}.')
    name, expression = enzymes.pop()
    if not name or not expression:
        raise ValueError('idXML has no usable digestion enzyme rule.')
    return {'name': name, 'regex': expression}


def main():
    """Write the OpenMS rule for the pepXML repair step."""
    args = parse_args()
    rule = extract_rule(args.input_path)
    with open(args.output_path, 'w', encoding='utf-8') as output:
        json.dump(rule, output)
    print(f'Extracted OpenMS enzyme rule: {rule["name"]} {rule["regex"]}')


if __name__ == '__main__':
    main()
