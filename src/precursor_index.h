#ifndef SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#define SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#include <vector>
#include "spectrum.h"

class precursor_index {

    // Contains all spectrum references, sorted first by charge, then by precursor mass
    std::vector<spectrum*> spectra;

public:

    explicit precursor_index(std::vector<spectrum*> &spectra);

    int get_size();
    int get_lower_bound(int charge, float mass);
    int get_upper_bound(int charge, float max_mass);
    float get_max_precursor_mass();
    spectrum *get_spectrum(int i);
private:

    int get_lower_bound_naive(int charge, float min_mass);
    int get_upper_bound_naive(int charge, float max_mass, int lower_bound=0);
};


#endif //SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
