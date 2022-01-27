#ifndef SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
#define SPECTRAL_SEARCH_ENGINE_SPECTRUM_H

#include <utility>
#include <vector>
#include <string>



class spectrum {
public:
    int id;
    std::string name;
    std::string peptide;
    //string species;
    float precursor_mass;
    int charge;

    int search_counter = 0;

    /*
     * Raw peaks
     * pos[i] corresponds to intensity[i]
     */

    std::vector<float> peak_positions;
    std::vector<float> intensities;

    /*
     * Binning and rescaling intensities
     */
    std::vector<float> bins; //vector over all bins (including zeros)
    std::vector<int> binned_peaks; //vector listing existing peaks binned
    std::vector<float> binned_intensities; //vector listing intensities corresponding to binned peaks
    int num_bins;

    //factor of intensity carried over to neighboring bins to account for mz-shifts
    float intensity_bin_spanning_factor = -0.5f; //set to -1.f to turn off (negative)
    bool remove_charge_reduced_precursor = false; //TODO uses spectrast magic function

    spectrum();

    bool bin_peaks(bool root_rescale=false, bool normalize=false);
    bool bin_peaks_sparse(bool root_rescale=false, bool normalize=false);

    bool denoise_mz_window(int topX, float window_size);

    bool root_scale_intensities();
    bool normalize_intensities();
    bool normalize_bins(float magnitude=-1.f);
    bool normalize_sparse_bins(float magnitude=-1.f);
    static int get_mz_bin(float mz);


    //Compare
    //friend bool operator<(const spectrum &one, const spectrum &other);
    bool operator<(const spectrum &other) const;
    bool operator<(std::pair<int,float> charge_mass_tuple) const;
    bool operator<=(std::pair<int,float> charge_mass_tuple) const;

private:
    bool add_intensity_to_bin(int bin, float intensity);
    bool spectrast_isNearPrecursor(double mz);
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
