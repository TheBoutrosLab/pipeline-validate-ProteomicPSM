"""Repair enzymatic-terminus annotations in OpenMS-exported pepXML.

OpenMS 3.5 writes num_tol_term=1 for every enzyme other than exactly
"Trypsin", while its PEP-derived probability is written only to the NTT=2
slot. This makes otherwise confident PSMs appear to have zero probability to
ProteinProphet. Keep the original PSM probability and calculate NTT for each
protein match from the original idXML enzyme rule and flanking residues.
"""

import argparse
from collections import Counter, namedtuple
import json
import re
import sys

from lxml import etree

# lxml exposes parser functions through a C extension.
# pylint: disable=c-extension-no-member


PEPXML_NS = 'http://regis-web.systemsbiology.net/pepXML'
NS = {'pep': PEPXML_NS}
N_TERMINAL = {'[', '-'}
C_TERMINAL = {']', '-'}
UNKNOWN = {'', 'X', '?'}
LOOKAROUND = re.compile(r'\(\?(?P<kind><=|<!|=|!)(?P<residues>\[[A-Z]+\]|[A-Z])\)')
EnzymeRule = namedtuple('EnzymeRule', 'name cleavage sense cut no_cut')


def parse_args():
    """Parse input and output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-path', required=True)
    parser.add_argument('--output-path', required=True)
    parser.add_argument('--enzyme-rule-path', required=True)
    return parser.parse_args()


def load_enzyme_rule(path):
    """Accept no cleavage or OpenMS rules with one-residue cleavage sites."""
    with open(path, encoding='utf-8') as rule_file:
        rule = json.load(rule_file)
    name, expression = rule['name'], rule['regex']
    if not name or not expression:
        raise ValueError('The idXML has no usable digestion enzyme rule.')

    # OpenMS represents no cleavage as an empty matching group. Protein
    # boundaries remain valid termini, but there are no internal cut sites.
    if name == 'no cleavage' and expression == '()':
        return EnzymeRule(name, None, 'C', '', '')

    assertions = list(LOOKAROUND.finditer(expression))
    if not assertions or ''.join(match.group() for match in assertions) != expression:
        raise ValueError(f'Unsupported OpenMS cleavage regex for {name}: {expression!r}.')
    positives = [match for match in assertions if match.group('kind') in ('<=', '=')]
    negatives = [match for match in assertions if match.group('kind') in ('<!', '!')]
    if len(positives) != 1 or len(negatives) > 1:
        raise ValueError(f'Unsupported OpenMS cleavage regex for {name}: {expression!r}.')
    sense = 'C' if positives[0].group('kind') == '<=' else 'N'
    if negatives and negatives[0].group('kind') != ('!' if sense == 'C' else '<!'):
        raise ValueError(f'Unsupported OpenMS cleavage regex for {name}: {expression!r}.')

    def residues(match):
        return match.group('residues').strip('[]')

    cut = residues(positives[0])
    no_cut = residues(negatives[0]) if negatives else ''
    return EnzymeRule(name, re.compile(expression), sense, cut, no_cut)


def update_pepxml_enzyme(run, enzyme):
    """Check the enzyme name and correct the converter's specificity fields."""
    pepxml_enzyme = run.find('pep:sample_enzyme', NS)
    if pepxml_enzyme is None:
        raise ValueError('pepXML run has no sample_enzyme.')
    rules = pepxml_enzyme.findall('pep:specificity', NS)
    if len(rules) != 1:
        raise ValueError('Expected exactly one enzyme specificity rule per run.')
    if pepxml_enzyme.get('name', '').casefold() != enzyme.name.casefold():
        raise ValueError(
            f'pepXML enzyme {pepxml_enzyme.get("name")!r} '
            f'differs from idXML enzyme {enzyme.name!r}.'
        )
    rules[0].set('sense', enzyme.sense)
    rules[0].set('cut', enzyme.cut)
    rules[0].set('no_cut', enzyme.no_cut)


def count_termini(sequence, before, after, cleavage):
    """Return NTT, or None when a flank is missing or unknown."""
    if not sequence or len(before) != 1 or len(after) != 1:
        return None
    if before in UNKNOWN or after in UNKNOWN:
        return None
    context = before + sequence + after
    n_term = before in N_TERMINAL or (
        cleavage is not None and cleavage.match(context, 1) is not None
    )
    c_term = after in C_TERMINAL or (
        cleavage is not None and cleavage.match(context, len(sequence) + 1) is not None
    )
    return int(n_term) + int(c_term)


def repair_hit(hit, cleavage, counts):
    """Repair the primary and alternative protein NTT values for one PSM."""
    sequence = hit.get('peptide', '')
    matches = [hit, *hit.findall('pep:alternative_protein', NS)]
    for match in matches:
        ntt = count_termini(
            sequence, match.get('peptide_prev_aa', ''),
            match.get('peptide_next_aa', ''), cleavage
        )
        if ntt is None:
            counts['unknown_flanks'] += 1
            continue
        if match.get('num_tol_term') != str(ntt):
            match.set('num_tol_term', str(ntt))
            counts['ntt_corrected'] += 1
        counts[f'ntt_{ntt}'] += 1

    results = hit.findall('pep:analysis_result/pep:peptideprophet_result', NS)
    if len(results) != 1:
        raise ValueError('Expected one peptideprophet_result per search_hit.')
    result = results[0]
    probability = result.get('probability')
    if probability is None or not 0 <= float(probability) <= 1:
        raise ValueError(f'Invalid PeptideProphet probability: {probability!r}.')
    # The PEP-derived score is not separately calibrated by NTT. Keep it
    # accessible for each protein match, whose termini may differ.
    result.set('all_ntt_prob', f'({probability},{probability},{probability})')
    counts['psms'] += 1


def repair_file(input_path, output_path, enzyme_rule_path):
    """Rewrite only enzyme-related pepXML attributes and preserve other data."""
    parser = etree.XMLParser(resolve_entities=False, no_network=True, huge_tree=True)
    tree = etree.parse(input_path, parser)
    root = tree.getroot()
    if root.tag != f'{{{PEPXML_NS}}}msms_pipeline_analysis':
        raise ValueError('Input is not a pepXML msms_pipeline_analysis document.')

    counts = Counter()
    runs = root.findall('pep:msms_run_summary', NS)
    if not runs:
        raise ValueError('pepXML has no msms_run_summary.')
    enzyme = load_enzyme_rule(enzyme_rule_path)
    for run in runs:
        update_pepxml_enzyme(run, enzyme)
        print(
            f'Repairing enzyme {enzyme.name}: sense={enzyme.sense}, '
            f'cut={enzyme.cut}, no_cut={enzyme.no_cut}'
        )
        for hit in run.findall('pep:spectrum_query/pep:search_result/pep:search_hit', NS):
            repair_hit(hit, enzyme.cleavage, counts)

    if not counts['psms']:
        raise ValueError('pepXML has no search_hit records.')
    tree.write(output_path, encoding='UTF-8', xml_declaration=True)
    print('Repair summary:', dict(sorted(counts.items())))
    if counts['unknown_flanks']:
        print(
            f'WARNING: Retained existing NTT for {counts["unknown_flanks"]} '
            'protein matches with missing or unknown flanks.',
            file=sys.stderr,
        )


def main():
    """Run the pepXML repair."""
    args = parse_args()
    repair_file(args.input_path, args.output_path, args.enzyme_rule_path)


if __name__ == '__main__':
    main()
