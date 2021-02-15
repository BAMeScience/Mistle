#ifndef SIMPLE_EXAMPLE_MSP_READER_H
#define SIMPLE_EXAMPLE_MSP_READER_H

#include "spectrum.h"

using namespace std;

enum msp_read_mode {
    DETAILED = 0,
    PEAKS,
    BINNED_PEAKS,
    SPARSE
};


class msp_reader {

static fstream infile;


public:

    static bool read_file(string &path, vector<spectrum*> &output_spectra, msp_read_mode read_mode=DETAILED);
    static bool read_file_precursors(string &path, vector<precursor *> &precursor_list);

    static bool read_spectra_from_positions(string &path, vector<precursor *> &precursor_list, vector<spectrum*> &output_spectra);
    static spectrum* read_spectrum_from_buffer(string buffer);


};


#endif //SIMPLE_EXAMPLE_MSP_READER_H
