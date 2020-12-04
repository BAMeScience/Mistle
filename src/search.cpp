#include <iostream>
#include <cmath>
#include "search.h"
#include "scores.h"

search::search() {

}

search::search(library *search_lib) : search_lib(search_lib) {

}

bool search::search_target_library(library *target_lib) {
    this->target_lib = target_lib;
    cout << "Begin searching target library" << endl;


    for (spectrum* query_spectrum : search_lib->spectrum_list) {
        float max_dot = 0.0;
        spectrum *best_candidate = nullptr;
        //Naive exhaustive search
        for (spectrum* candidate_spectrum : target_lib->spectrum_list) {
            if (is_candidate_suitable(candidate_spectrum, query_spectrum)) {
                float dot = scores::dot_product(query_spectrum->bins, candidate_spectrum->bins);
                if (dot >= max_dot) {//todo What if equal
                    best_candidate = candidate_spectrum;
                    max_dot = dot;
                }
            }
        }
        match best_match(query_spectrum, best_candidate, max_dot, 1);
        search_results.push_back(best_match);
    }

    return true;
}

bool search::is_candidate_suitable(spectrum *candidate_spectrum, spectrum *query_spectrum) {
    bool has_equal_charge = candidate_spectrum->charge == query_spectrum->charge;
    bool is_in_mass_range = abs(candidate_spectrum->precursor_mass - query_spectrum->precursor_mass) < mz_tolerance;
    return has_equal_charge && is_in_mass_range;
}

