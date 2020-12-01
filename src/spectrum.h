#ifndef SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
#define SPECTRAL_SEARCH_ENGINE_SPECTRUM_H

#include <utility>
#include <vector>
#include <string>

using namespace std;

class spectrum {
public:
    string peptide;
    float precursor_mass;
    int charge;

    // peak entries correspond to intensities 1 to 1 at each position.
    vector<float> peaks;
    vector<float> intensities;

    vector<float> bins;

    spectrum();

private:
    bool bin_peaks();
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
