#ifndef SIMPLE_EXAMPLE_SEARCH_H
#define SIMPLE_EXAMPLE_SEARCH_H
#include "library.h"
#include "match.h"

using namespace std;

class search {
    library *search_lib;
    library *target_lib;

    vector<match> search_results;

public:
    search();
    search(library *search_lib);

    bool search_target_library(library *target_lib);
};


#endif //SIMPLE_EXAMPLE_SEARCH_H
