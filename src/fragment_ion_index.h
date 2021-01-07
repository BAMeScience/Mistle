#ifndef SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#define SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#include <vector>
#include "precursor_index.h"

using namespace std;

struct fragment {
    int parent_id;
    float intensity;

    fragment(int parent_id, float intensity) : parent_id(parent_id), intensity(intensity) {};
};


typedef vector<fragment> fragment_bin;


class fragment_ion_index {
public:
    vector<fragment_bin> fragment_bins;

    explicit fragment_ion_index(precursor_index *parent_index);
};


#endif //SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
