#ifndef SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
#define SPECTRAL_SEARCH_ENGINE_SPECTRUM_H

#include <utility>
#include <vector>
#include <string>



struct precursor {
    int id = 1000;
    int rank = 1000;
    float mass;
    int charge;

    unsigned long offset_begin;
    unsigned long offset_end;
    std::string name;
    //string peptide;
    precursor() {};
    precursor(int id, float mass, int charge, std::string name, std::string peptide="") : id(id), mass(mass), charge(charge) {};

    bool operator<(const precursor &other) const {
        return charge < other.charge || (charge == other.charge && mass < other.mass);
    };

};

class spectrum {
public:
    int id;
    std::string name;
    std::string peptide;
    //string species;
    float precursor_mass;
    int charge;


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
    bool operator<(std::pair<int,float> charge_mass_tuple) const;
    bool operator<=(std::pair<int,float> charge_mass_tuple) const;

private:
    static float rescale_intensity(float intensity);
    bool spectrast_isNearPrecursor(double mz);
};


#endif //SPECTRAL_SEARCH_ENGINE_SPECTRUM_H
