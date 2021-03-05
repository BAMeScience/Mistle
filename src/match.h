#ifndef SIMPLE_EXAMPLE_MATCH_H
#define SIMPLE_EXAMPLE_MATCH_H

#include "spectrum.h"


class match {
public:
    spectrum *query_spectrum;
    spectrum *matched_spectrum;

    unsigned int query_id;
    unsigned int target_id;

    float dot_product;
    float mass_difference;
    int hit_rank;
    int query_index; //todo implement

    int num_matched_peaks;

    match();
    match(spectrum *search_spectrum, spectrum *matched_spectrum);
    match(spectrum *search_spectrum, spectrum *matched_spectrum, float dot_product, int hit_rank);
    match(unsigned int query_id, unsigned int target_id, float dot_product, float mass_difference, int hit_rank);

};


#endif //SIMPLE_EXAMPLE_MATCH_H
