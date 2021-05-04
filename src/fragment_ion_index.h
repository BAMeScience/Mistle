#ifndef SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#define SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
#include <vector>
#include <immintrin.h>
#include "precursor_index.h"


struct fragment {
    unsigned int parent_id;
    float intensity;
    float mz;

    // If multiple peaks are contributing to the fragment - keep track of composition
    std::vector<std::pair<float,float>> peak_composition; //<MZ, INTENSITY>

    fragment(unsigned int parent_id, float intensity) : parent_id(parent_id), intensity(intensity) {};
    fragment(unsigned int parent_id, float intensity, float mz) : parent_id(parent_id), intensity(intensity), mz(mz) {};
};


typedef std::vector<fragment> fragment_bin;



struct __attribute__ ((aligned (32))) fragment_binn {
    __attribute__ ((aligned (32))) std::vector<float> intensities;
    __attribute__ ((aligned (32))) std::vector<unsigned> parent_ids;
#if USE_AVX_512
    __attribute__ ((aligned (32))) std::vector<__m512> _intensities;
    __attribute__ ((aligned (32))) std::vector<__m512i> _parent_ids;

#elif USE_AVX_2
    __attribute__ ((aligned (32))) std::vector<__m256> _intensities;
    __attribute__ ((aligned (32))) std::vector<__m256i> _parent_ids;
    __attribute__ ((aligned (32))) std::vector<__m256i> _parent_ranks;
#endif
};

class fragment_ion_index {
public:

    std::shared_ptr<precursor_index> precursor_idx;
    std::string file_path;
    std::vector<fragment_bin> fragment_bins;
    __attribute__ ((aligned (32))) std::vector<fragment_binn> frag_bins;


    fragment_ion_index();
    explicit fragment_ion_index(precursor_index *parent_index);
    explicit fragment_ion_index(std::string path);

    bool sort_index(std::unique_ptr<precursor_index>& parent_index);


    bool prepare_axv_access();
    bool load_index_from_file(const std::string& path);
    bool load_index_from_binary_file(const std::string& path);
    bool load_preliminary_index_from_binary_file(const std::string& path);
    bool save_index_to_file(const std::string& path);
    bool save_index_to_binary_file(const std::string& path);

};


#endif //SIMPLE_EXAMPLE_FRAGMENT_ION_INDEX_H
