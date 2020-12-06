#include "match.h"

match::match() {

}

match::match(spectrum *search_spectrum, spectrum *matched_spectrum, float dot_product, int hit_rank) : query_spectrum(search_spectrum), matched_spectrum(matched_spectrum), dot_product(dot_product), hit_rank(hit_rank) {

}
