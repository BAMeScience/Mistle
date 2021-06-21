#ifndef SIMPLE_EXAMPLE_MSP_READER_H
#define SIMPLE_EXAMPLE_MSP_READER_H

#include "spectrum.h"
#include "precursor_index.h"
#include <memory>


class msp_reader {

static std::fstream infile;


public:

    static bool read_file(std::string &path, std::vector<std::shared_ptr<spectrum>> &output_spectra);
    static bool continue_read_file(std::vector<std::shared_ptr<spectrum>> &output_spectra);
    static bool read_file_precursors(std::string &path, std::vector<precursor *> &precursor_list);
    static bool read_file_precursors_efficient(std::string &path, std::vector<precursor *> &precursor_list);

    static bool read_spectra_from_positions(std::string &path, std::vector<precursor *> &precursor_list, std::vector<spectrum*> &output_spectra);
    static bool read_next_entry_into_buffer(std::ifstream &f, std::string &buffer);
    static std::shared_ptr<spectrum> read_spectrum_from_buffer(const std::string& buffer);


};


#endif //SIMPLE_EXAMPLE_MSP_READER_H
