#ifndef SIMPLE_EXAMPLE_SPECTRAL_SEARCH_H
#define SIMPLE_EXAMPLE_SPECTRAL_SEARCH_H
#include "library.h"
#include "match.h"

using namespace std;

class spectral_search {
    library *search_lib;
    library *target_lib;

    vector<match> search_results;
    float mz_tolerance=3.0;

public:
    spectral_search();
    explicit spectral_search(library *search_lib);

    bool search_target_library(library *target_lib);
    vector<match> search_specific_by_name(string query_name, library *target_lib, int num_matches);
    vector<match> get_results();
    bool save_results_to_file(string path, string delimiter="\t");

private:
    bool is_candidate_suitable(spectrum *candidate_spectrum, spectrum *query_spectrum); //Checking if charge and mass conditions are fullfilled to warrant a closer look
};


#endif //SIMPLE_EXAMPLE_SPECTRAL_SEARCH_H
