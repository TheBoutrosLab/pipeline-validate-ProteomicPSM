""" This script reads in an idXML file and fix the `spectra_file` UserParam
field of the ID run. This is just to resolve the issue the ConsensusID only
runs with ProteinIdentification with exactly one mzML file. """
import argparse
from pyopenms import IdXMLFile, PeptideIdentificationList

def parse_args() -> argparse.Namespace:
    """ parse args """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--input-path',
        type=str,
        help='Path to the input idXML file.',
        required=True
    )
    parser.add_argument(
        '--output-path',
        type=str,
        help='Path to the output idXML file.',
        required=True
    )
    parser.add_argument(
        '--mzML',
        type=str,
        help='Path to the mzML file.',
        required=True
    )
    return parser.parse_args()

def main():
    """ File spectra file """
    args = parse_args()

    prot_idents = []
    pep_idents = PeptideIdentificationList()
    IdXMLFile().load(args.input_path, prot_idents, pep_idents)

    if len(prot_idents) > 1:
        raise ValueError('The idXML contains more than 1 run.')

    prot_ident = prot_idents[0]
    mzml_file = args.mzML.encode('utf-8')
    prot_ident.setPrimaryMSRunPath([mzml_file], False)

    IdXMLFile().store(args.output_path, [prot_ident], pep_idents)


if __name__ == "__main__":
    main()
