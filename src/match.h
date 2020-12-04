#ifndef SIMPLE_EXAMPLE_MATCH_H
#define SIMPLE_EXAMPLE_MATCH_H

#include "spectrum.h"

using namespace std;

class match {
public:
    spectrum *search_spectrum;
    spectrum *matched_spectrum;

    float dot_product;
    int hit_rank;

    match();
    match(spectrum *search_spectrum, spectrum *matched_spectrum, float dot_product, int hit_rank);

};


#endif //SIMPLE_EXAMPLE_MATCH_H
