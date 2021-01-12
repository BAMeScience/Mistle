#ifndef SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
#define SPECTRAL_SEARCH_ENGINE_SPECTRUM_H

#include <utility>
#include <vector>
#include <string>


using namespace std;

class spectrum {
public:
    int id;
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
    vector<float> bins; //vector over all bins (including zeros)
    vector<int> binned_peaks; //vector listing existing peaks binned
    vector<float> binned_intensities; //vector listing intensities corresponding to binned peaks


    //factor of intensity carried over to neighboring bins to account for mz-shifts
    float intensity_bin_spanning_factor = -0.5f; //set to -1.f to turn off (negative)
    bool remove_charge_reduced_precursor = true; //TODO uses spectrast magic function

    spectrum();

    bool bin_peaks(bool root_rescale=false, bool normalize=false);
    bool bin_peaks_sparse(bool root_rescale=false, bool normalize=false);
    bool normalize_bins(float magnitude=-1.f);
    bool normalize_sparse_bins(float magnitude=-1.f);

    //Compare
    //friend bool operator<(const spectrum &one, const spectrum &other);
    bool operator<(const spectrum &other) const;
    bool operator<(pair<int,float> charge_mass_tuple) const;
    bool operator<=(pair<int,float> charge_mass_tuple) const;

private:
    static float rescale_intensity(float intensity);
    bool spectrast_isNearPrecursor(double mz);
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
