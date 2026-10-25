"""Check pepXML terminal counts for an enzyme with no internal cuts."""

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from xml.etree import ElementTree


ROOT = Path(__file__).resolve().parents[1]
PEPXML_NS = 'http://regis-web.systemsbiology.net/pepXML'
NS = {'pep': PEPXML_NS}


class PepXMLTerminiTest(unittest.TestCase):
    """Exercise the repair with protein-boundary and internal peptide hits."""

    @classmethod
    def setUpClass(cls):
        """Load the repair script without requiring it to be a package."""
        spec = importlib.util.spec_from_file_location(
            'FixPepXMLTermini_test', ROOT / 'src/FixPepXMLTermini.py'
        )
        cls.repair = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.repair)

    def test_no_cleavage_counts_only_protein_boundaries(self):
        """Keep PSM probabilities while assigning NTT 2, 1, and 0."""
        queries = []
        for scan, before, after in ((1, '-', '-'), (2, '-', 'A'), (3, 'A', 'A')):
            queries.append(
                f'<spectrum_query spectrum="sample.{scan}.{scan}.2">'
                '<search_result><search_hit peptide="PEPTIDE" '
                f'peptide_prev_aa="{before}" peptide_next_aa="{after}" '
                'num_tol_term="1">'
                '<analysis_result analysis="peptideprophet">'
                '<peptideprophet_result probability="0.95" '
                'all_ntt_prob="(0,0,0.95)"/>'
                '</analysis_result></search_hit></search_result></spectrum_query>'
            )
        source = (
            f'<msms_pipeline_analysis xmlns="{PEPXML_NS}">'
            '<msms_run_summary><sample_enzyme name="no cleavage">'
            '<specificity cut="" sense="C"/></sample_enzyme>'
            f'{"".join(queries)}</msms_run_summary></msms_pipeline_analysis>'
        )
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            input_path = directory / 'input.pepXML'
            output_path = directory / 'output.pepXML'
            rule_path = directory / 'rule.json'
            input_path.write_text(source, encoding='utf-8')
            rule_path.write_text(
                json.dumps({'name': 'no cleavage', 'regex': '()'}), encoding='utf-8'
            )
            self.repair.repair_file(str(input_path), str(output_path), str(rule_path))
            root = ElementTree.parse(output_path).getroot()

        hits = root.findall('.//pep:search_hit', NS)
        self.assertEqual([hit.get('num_tol_term') for hit in hits], ['2', '1', '0'])
        for hit in hits:
            result = hit.find('pep:analysis_result/pep:peptideprophet_result', NS)
            self.assertEqual(result.get('probability'), '0.95')
            self.assertEqual(result.get('all_ntt_prob'), '(0.95,0.95,0.95)')
        specificity = root.find('.//pep:sample_enzyme/pep:specificity', NS)
        self.assertEqual(specificity.get('cut'), '')


if __name__ == '__main__':
    unittest.main()
