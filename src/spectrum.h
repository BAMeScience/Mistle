#ifndef SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
#define SPECTRAL_SEARCH_ENGINE_SPECTRUM_H

#include <utility>
#include <vector>
#include <string>


using namespace std;

class spectrum {
public:
    string name;
    string peptide;
    float precursor_mass;
    int charge;

    // peak positions correspond to intensities 1 to 1 for each vector entry.
    vector<float> peak_positions;
    vector<float> intensities;

    vector<float> bins;
    //factor of intensity carried over to neighboring bins to account for mz-shifts
    float intensity_bin_spanning_factor = 0.5f; //set to -1.f to turn off

    spectrum();

    bool bin_peaks(bool root_rescale=false, bool normalize=false);
    bool normalize_bins(float magnitude=-1);

private:
    static float rescale_intensity(float intensity);
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
