#ifndef SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#define SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
#include <vector>
#include "spectrum.h"

using namespace std;

class precursor_index {
public:
    vector<spectrum*> spectra;

    explicit precursor_index(vector<spectrum*> &spectra);

    int get_lower_bound(float min_mass);
    int get_upper_bound(float max_mass);
};


#endif //SIMPLE_EXAMPLE_PRECURSOR_INDEX_H
