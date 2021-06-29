#ifndef SIMPLE_EXAMPLE_MGF_READER_H
#define SIMPLE_EXAMPLE_MGF_READER_H

#include <memory>
#include "spectrum.h"


class mgf_reader {
public:
    mgf_reader();

    static bool read_file(std::string path, std::vector<std::shared_ptr<spectrum>> &output_spectra);
    static bool read_file_batch(std::fstream &infile, std::vector<std::shared_ptr<spectrum>> &output_spectra, int batch_size);
};


#endif //SIMPLE_EXAMPLE_MGF_READER_H
