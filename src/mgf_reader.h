#ifndef SIMPLE_EXAMPLE_MGF_READER_H
#define SIMPLE_EXAMPLE_MGF_READER_H
#include "spectrum.h"


class mgf_reader {
public:
    mgf_reader();

    static bool read_file(std::string path, std::vector<spectrum *> &output_spectra);
};


#endif //SIMPLE_EXAMPLE_MGF_READER_H
