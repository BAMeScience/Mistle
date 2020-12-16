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


    /*
     * Raw peaks
     * pos[i] corresponds to intensity[i]
     */

    vector<float> peak_positions;
    vector<float> intensities;

    /*
     * Binning and rescaling intensities
     */
    vector<float> bins;
    //factor of intensity carried over to neighboring bins to account for mz-shifts
    float intensity_bin_spanning_factor = -0.5f; //set to -1.f to turn off
    bool remove_charge_reduced_precursor = true; //TODO uses spectrast magic function

    spectrum();

    bool bin_peaks(bool root_rescale=false, bool normalize=false);
    bool normalize_bins(float magnitude=-1.f);

private:
    static float rescale_intensity(float intensity);
    bool spectrast_isNearPrecursor(double mz);
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
