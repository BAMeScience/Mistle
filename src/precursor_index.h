#ifndef SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#define SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#include <vector>
#include "spectrum.h"

using namespace std;

class precursor_index {

    vector<spectrum*> spectra; //todo one for each charge

public:

    explicit precursor_index(vector<spectrum*> &spectra);

    int get_size();
    int get_lower_bound(float min_mass);
    int get_upper_bound(float max_mass);
    float get_max_precursor_mass();
    spectrum *get_spectrum(int i);
};


#endif //SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
