#ifndef SIMPLE_EXAMPLE_SPECTRAL_SEARCH_H
#define SIMPLE_EXAMPLE_SPECTRAL_SEARCH_H
#include "library.h"
#include "match.h"
#include "fragment_ion_index.h"


class spectral_search {
    library *query_lib;
    library *target_lib;

    std::vector<match> search_results;
    float mz_tolerance=3.0;

public:
    spectral_search();
    explicit spectral_search(library *query_lib);
    spectral_search(library *query_lib, library *target_lib);

    /*
     * Searching query library against a target library
     */
    bool search_target_library();
    bool search_target_library(library *target_lib);
    bool search_fragment_ion_index();

    std::vector<match> get_results();
    bool save_results_to_file(std::string path, std::string delimiter="\t");
    bool read_results_from_file(std::string path, char delimiter='\t', bool read_dot=false, bool has_header=true);

    /*
     * (Re-)scoring of matches
     */
    bool rescore_matches();

private:
    bool is_candidate_suitable(spectrum *candidate_spectrum, spectrum *query_spectrum); //Checking if charge and mz conditions are fullfilled to warrant a closer look
};


#endif //SIMPLE_EXAMPLE_SPECTRAL_SEARCH_H
