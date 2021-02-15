#ifndef SIMPLE_EXAMPLE_LIBRARY_H
#define SIMPLE_EXAMPLE_LIBRARY_H
#include <vector>
#include "spectrum.h"
#include "precursor_index.h"
#include "fragment_ion_index.h"


class library {
public:
    library();
    library(std::string &path);
    library(std::vector<spectrum*> &spectra);
    ~library();

    bool load_spectra_from_file(std::string path);
    bool load_library_from_directory(std::string &path);

    bool build_library_index();

    std::vector<spectrum*> spectrum_list;
    precursor_index* precursor_idx;
    fragment_ion_index* fragment_ion_idx;
    bool is_indexed = false;

private:

};


#endif //SIMPLE_EXAMPLE_LIBRARY_H
