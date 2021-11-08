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

    float mass_difference;

    // Default scores
    float similarity;
    float bias;
    float dot_product;

    //Additional score
    float delta_dot;
    float delta_similarity;
    int peak_count_query;
    int peak_count_target;


    //Advanced scores
    float sim2;
    float spectraST_score;
    float spectraST_score_dot;
    double x_hunter_score; //all peaks, not just top 20
    double x_hunter_score_dot; //all peaks, not just top 20


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
