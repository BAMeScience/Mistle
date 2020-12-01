#ifndef SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
#define SPECTRAL_SEARCH_ENGINE_SPECTRUM_H

#include <utility>
#include <vector>
#include <string>

using namespace std;

class spectrum {
    vector<pair<float, float>> peak_intensities;
    vector<float> unit_bins;

    float mass;
    string peptide;

public:
    spectrum();

private:
    bool bin_peaks();
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
