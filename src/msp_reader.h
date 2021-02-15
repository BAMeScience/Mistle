#ifndef SIMPLE_EXAMPLE_MSP_READER_H
#define SIMPLE_EXAMPLE_MSP_READER_H

#include "spectrum.h"


enum msp_read_mode {
    DETAILED = 0,
    PEAKS,
    BINNED_PEAKS,
    SPARSE
};


class msp_reader {

static std::fstream infile;


public:

    static bool read_file(std::string &path, std::vector<spectrum*> &output_spectra, msp_read_mode read_mode=DETAILED);
    static bool read_file_precursors(std::string &path, std::vector<precursor *> &precursor_list);
    static bool read_file_precursors_efficient(std::string &path, std::vector<precursor *> &precursor_list);
    static bool read_file_precursors_efficient2(std::string &path, std::vector<precursor *> &precursor_list);

    static bool read_spectra_from_positions(std::string &path, std::vector<precursor *> &precursor_list, std::vector<spectrum*> &output_spectra);
    static spectrum* read_spectrum_from_buffer(std::string buffer);


};


#endif //SIMPLE_EXAMPLE_MSP_READER_H
