#include "match.h"

match::match() {

}

match::match(spectrum *search_spectrum, spectrum *matched_spectrum, float dot_product, int hit_rank) : query_spectrum(search_spectrum), matched_spectrum(matched_spectrum), dot_product(dot_product), mass_difference(mass_difference), hit_rank(hit_rank) {
    mass_difference = matched_spectrum->precursor_mass - search_spectrum->precursor_mass;
}
