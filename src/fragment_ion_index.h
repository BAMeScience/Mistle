#ifndef SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#define SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#include <vector>
#include "precursor_index.h"


struct fragment {
    int parent_id;
    float intensity;

    fragment(int parent_id, float intensity) : parent_id(parent_id), intensity(intensity) {};
};


typedef std::vector<fragment> fragment_bin;


class fragment_ion_index {
public:
    std::string file_path;
    std::vector<fragment_bin> fragment_bins;

    explicit fragment_ion_index(precursor_index *parent_index);
    explicit fragment_ion_index(std::string path);

    bool sort_index(precursor_index &parent_index);



    bool load_index_from_file(const std::string& path);
};


#endif //SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
