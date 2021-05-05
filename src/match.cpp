#include "match.h"
#include <memory>
#include <iostream>

using namespace std;

match::match() {

}

match::match(std::shared_ptr<spectrum> search_spectrum, std::shared_ptr<spectrum> matched_spectrum) : query_spectrum(search_spectrum), matched_spectrum(matched_spectrum) {
    mass_difference = matched_spectrum->precursor_mass - search_spectrum->precursor_mass;
}

match::match(std::shared_ptr<spectrum> search_spectrum, std::shared_ptr<spectrum> matched_spectrum, float dot_product, int hit_rank) : query_spectrum(search_spectrum), matched_spectrum(matched_spectrum), dot_product(dot_product), hit_rank(hit_rank) {
    mass_difference = matched_spectrum->precursor_mass - search_spectrum->precursor_mass;
}

match::match(unsigned int query_id, unsigned int target_id, float dot_product, float mass_difference, int hit_rank) : query_id(query_id), target_id(target_id), dot_product(dot_product), mass_difference(mass_difference), hit_rank(hit_rank) {

}

match::match(unsigned int query_id, unsigned int target_id, float similarity_score, float dot_product,
             float mass_difference, int hit_rank) : query_id(query_id), target_id(target_id), similarity_score(similarity_score), dot_product(dot_product), mass_difference(mass_difference), hit_rank(hit_rank) {

}
