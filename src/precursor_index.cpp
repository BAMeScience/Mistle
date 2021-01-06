#include "precursor_index.h"

precursor_index::precursor_index(vector<spectrum *> &spectra) : spectra(spectra) {
    sort(this->spectra.begin(), this->spectra.end()); //TODO add comparator

    for (int i = 0; i < spectra.size(); ++i) {
        spectra[i]->id = i;
    }

}

//TODO replace with binary search
int precursor_index::get_lower_bound(float min_mass) {
    for (int i = 0; i < spectra.size(); ++i) {
        if (spectra[i]->precursor_mass >= min_mass) {
            return i;
        }
    }
    return -1;
}

//TODO replace with binary search
int precursor_index::get_upper_bound(float max_mass) {
    for (int i = 0; i < spectra.size(); ++i) {
        if (spectra[i]->precursor_mass > max_mass) {
            return i-1;
        }
    }
}
