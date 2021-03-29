#ifndef SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#define SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#include <vector>
#include "precursor_index.h"


struct fragment {
    unsigned int parent_id;
    float intensity;

    fragment(unsigned int parent_id, float intensity) : parent_id(parent_id), intensity(intensity) {};
};


typedef std::vector<fragment> fragment_bin;



struct __attribute__ ((aligned (32))) fragment_binn {
    __attribute__ ((aligned (32))) std::vector<float> intensities;
};

class fragment_ion_index {
public:
    std::string file_path;
    std::vector<fragment_bin> fragment_bins;
    __attribute__ ((aligned (32))) std::vector<fragment_binn> frag_bins;


    fragment_ion_index();
    explicit fragment_ion_index(precursor_index *parent_index);
    explicit fragment_ion_index(std::string path);

    bool sort_index(std::unique_ptr<precursor_index>& parent_index);


    bool update_intensities();
    bool load_index_from_file(const std::string& path);
    bool load_index_from_binary_file(const std::string& path);
    bool save_index_to_file(const std::string& path);
    bool save_index_to_binary_file(const std::string& path);
};


#endif //SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
