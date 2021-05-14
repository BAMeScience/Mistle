#ifndef SIMPLE_EXAMPLE_MATCH_H
#define SIMPLE_EXAMPLE_MATCH_H

#include "spectrum.h"
#include <memory>

class match {
public:
    std::shared_ptr<spectrum> query_spectrum;
    std::shared_ptr<spectrum> matched_spectrum;

    unsigned int query_id;
    unsigned int target_id;

    float similarity_score;
    float dot_product;
    float mass_difference;
    int hit_rank;
    int query_index; //todo implement

    int num_matched_peaks;

    match();
    match(unsigned int query_id, unsigned int target_id);

    /*
     * Deprecated constructors
     */
    match(std::shared_ptr<spectrum> search_spectrum, std::shared_ptr<spectrum> matched_spectrum);
    match(std::shared_ptr<spectrum> search_spectrum, std::shared_ptr<spectrum> matched_spectrum, float dot_product, int hit_rank);
    match(unsigned int query_id, unsigned int target_id, float dot_product, float mass_difference, int hit_rank);
    match(unsigned int query_id, unsigned int target_id, float similarity_score, float dot_product, float mass_difference, int hit_rank);

};


#endif //SIMPLE_EXAMPLE_MATCH_H
