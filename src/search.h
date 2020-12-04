#ifndef SIMPLE_EXAMPLE_SEARCH_H
#define SIMPLE_EXAMPLE_SEARCH_H
#include "library.h"
#include "match.h"

using namespace std;

class search {
    library *search_lib;
    library *target_lib;

    vector<match> search_results;
    float mz_tolerance=3.0;

public:
    search();
    search(library *search_lib);

    bool search_target_library(library *target_lib);

private:
    bool is_candidate_suitable(spectrum *candidate_spectrum, spectrum *query_spectrum); //Checking if charge and mass conditions are fullfilled to warrant a closer look
};


#endif //SIMPLE_EXAMPLE_SEARCH_H
