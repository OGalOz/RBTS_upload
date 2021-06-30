/*
A KBase module: RBTS_upload
*/

module RBTS_upload {
    typedef structure {
        string report_name;
        string report_ref;
    } ReportResults;

    /*
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_RBTS_upload(mapping<string,UnspecifiedObject> params) returns (ReportResults output) authentication required;

};
