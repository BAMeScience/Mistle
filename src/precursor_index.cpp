#include <iostream>
#include <utility>
#include "precursor_index.h"

precursor_index::precursor_index(vector<spectrum *> &spectrum_list) : spectra(spectrum_list) {
    sort(spectra.begin(), spectra.end(), [](const spectrum *a, const spectrum *b) {
        return *a < *b;
    });

    for (int i = 0; i < spectra.size(); ++i) {
        spectra[i]->id = i;
    }
}

//TODO replace with binary search
int precursor_index::get_lower_bound(int charge, float min_mass) {
    for (int i = 0; i < spectra.size(); ++i) {
        if (spectra[i]->charge == charge && spectra[i]->precursor_mass >= min_mass) {
            return i; //Find first precursor with matching charge and greater mass
        }
    }
    return -1;
}

//TODO replace with binary search
int precursor_index::get_upper_bound(int charge, float max_mass, int lower_bound) {
    for (int i = lower_bound; i < spectra.size(); ++i) {
        if (spectra[i]->charge > charge || (spectra[i]->charge == charge && spectra[i]->precursor_mass > max_mass)) {
            return i-1; //Find last elem, that has not higher charge and does not exceed max mass with matching charge
        }
    }
    return spectra.size() - 1;
}

int precursor_index::get_size() {
    return spectra.size();
}

spectrum *precursor_index::get_spectrum(int i) {
    return spectra[i];
}

float precursor_index::get_max_precursor_mass() {
    return spectra.back()->precursor_mass;
}
