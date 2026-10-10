include { convert_idXML2pepXML } from "${moduleDir}/convert_idXML2pepXML"
include { run_PeptideProphet } from "${moduleDir}/run_PeptideProphet"
include { run_iProphet } from "${moduleDir}/run_iProphet"
include { run_ProteinProphet } from "${moduleDir}/run_ProteinProphet"
include { filter_PSM_philosopher } from "${moduleDir}/filter_PSM_philosopher"
include { convert_pepXML2idXML } from "${moduleDir}/convert_pepXML2idXML"
include { filter_PSM } from "${moduleDir}/filter_PSM"

workflow validate_PSM_philosopher {
    take:
        ich
    main:
        convert_idXML2pepXML(ich)
        run_PeptideProphet(convert_idXML2pepXML.out[0], file(params.decoy_database))
        pepXMLs = run_PeptideProphet.out[0].collect()
        run_iProphet(pepXMLs, file(params.decoy_database))
        // run_ProteinProphet(run_iProphet.out[0], file(params.decoy_database))
        // protXML = run_ProteinProphet.out[0]
        convert_pepXML2idXML(run_iProphet.out[0])
        filter_PSM(convert_pepXML2idXML.out[0])
}
