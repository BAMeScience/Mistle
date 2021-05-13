#include <iostream>
#include <cmath>
#include <numeric>
#include "spectrum.h"
#include "DefineConstants.h"
#include "settings.h"

using namespace std;


spectrum::spectrum() {

}


bool spectrum::bin_peaks(bool root_rescale, bool normalize) {

    //TODO: OBSOLETE use sparse binning instead

    num_bins = int((BIN_MAX_MZ - BIN_MIN_MZ) / settings::bin_size) + 1; // TODO maybe smart

    bins = vector<float>(num_bins, 0.f);

    for (int i = 0; i < peak_positions.size(); ++i) {
        int bin = get_mz_bin(peak_positions[i]);
        if (bin < 0 || bin > bins.size() - 1) { // TODO spectraST light ion, cut_off. bin < 180 ??
            //TODO what would spectrast do?
            //cout << "Warning peak out of bin range :: discarding intensity" << endl;
            bin = 0;

            continue;
        }
        /*if (abs(float(bin) - precursor_mass) < 10) //TODO why is this a thing?
            continue;*/
        if (remove_charge_reduced_precursor && spectrast_isNearPrecursor(peak_positions[i])) {
            continue;
        }

        float intensity = intensities[i];
        if (root_rescale)
            intensity = sqrt(intensity);

        bins[bin] = sqrt(bins[bin] * bins[bin] + intensity + intensity); //sqrt-accumulate if multiple peaks fall into the sam bin
        if (intensity_bin_spanning_factor > 0.f) {
            float neighbor_intensity = intensity * intensity_bin_spanning_factor;
            if (bin > 0) {
                bins[bin-1] = sqrt(intensity * intensity + neighbor_intensity * neighbor_intensity);
            }
            if (bin < bins.size()) {
                bins[bin+1] = sqrt(intensity * intensity + neighbor_intensity * neighbor_intensity);
            }
        }
    }

    if (normalize) {
        /*float min_cut = 0.01f; //TODO undirty or delete this
        for (float &i : bins) {
            if (i<min_cut) {
                i = 0.f;
            }
        }*/
        return normalize_bins();
    }
    return true;
}


bool spectrum::spectrast_isNearPrecursor(double mz) {

    if (precursor_mass < 0.0001) return (false); // m_parentMz not set, can't determine

    /*if (m_fragType != "ETD" && m_fragType != "ETD-SA" && m_fragType != "ETD-HR") {
        if (mz > m_parentMz - 60.0 && mz < m_parentMz + 20.0) {
            return (true);
        }
    } else {*/

    int lowCharge = charge; //TODO this is the tool I am dealing with
    int highCharge = charge;
    //    if (m_parentCharge == 0) {
    // remove all possible charge-reduced precursors for precursor charge up to 6.
    lowCharge = 3;
    highCharge = 6;
    //}

    for (int parentCharge = lowCharge; parentCharge <= highCharge; parentCharge++) {
        double parentMass = precursor_mass * parentCharge;
        for (int c = parentCharge; c >= 1; c--) {
            if (mz >= (parentMass - 20.0) / (double)c && mz <= (parentMass + 6.0) / (double)c) {
                return (true);
            }
        }
    }
    return (false);
}


bool spectrum::normalize_bins(float magnitude) {
    if (magnitude < 0.f) { // Default, magnitude not specified
        magnitude = 0.f;
        for (float &i : bins) {
            magnitude += i*i;
        }
        magnitude = sqrt(magnitude);
    }
    for (float &i : bins) {
        i = i / magnitude;
    }

    return true;
}

bool spectrum::bin_peaks_sparse(bool root_rescale, bool normalize) {

    binned_peaks.clear();
    binned_intensities.clear();
    num_bins = int((BIN_MAX_MZ - BIN_MIN_MZ) / settings::bin_size) + 1;

    for (int i = 0; i < peak_positions.size(); ++i) {

        if (peak_positions[i] < BIN_MIN_MZ || peak_positions[i] > BIN_MAX_MZ) {
            continue;
        }

        if (remove_charge_reduced_precursor && spectrast_isNearPrecursor(peak_positions[i])) {
            continue;
        }
        //Determine bin
        int bin = get_mz_bin(peak_positions[i]);

        //Retrieve and scale intensity
        float intensity = intensities[i];
        if (root_rescale)
            intensity = sqrt(intensity);

        //Update existing bin (if possible)
        add_intensity_to_bin(bin, intensity);
        for (int j = 1; j <= settings::neighbors; ++j) { //Adding intensity fraction to neighboring bins
            add_intensity_to_bin(bin + j, intensity * pow(settings::neighbors_intensity_factor, j));
            add_intensity_to_bin(bin - j, intensity * pow(settings::neighbors_intensity_factor, j));
        }

    }

    //Normalize sparse bins
    if (normalize) {
        return normalize_sparse_bins();
    }
    return true;
}

bool spectrum::normalize_sparse_bins(float magnitude) {
    if (magnitude < 0.f) { // Default, magnitude not specified
        magnitude = 0.f;
        for (float &i : binned_intensities) {
            magnitude += i*i;
        }
        magnitude = sqrt(magnitude);
    }
    for (float &i : binned_intensities) {
        i = i / magnitude;
    }

    return true;
}

bool spectrum::operator<(const spectrum &other) const {
    return charge < other.charge || (charge == other.charge && precursor_mass < other.precursor_mass);
}

bool spectrum::operator<(pair<int, float> charge_mass_tuple) const {
    return charge < charge_mass_tuple.first || (charge == charge_mass_tuple.first && precursor_mass < charge_mass_tuple.second);
}

bool spectrum::operator<=(pair<int, float> charge_mass_tuple) const {
    return charge < charge_mass_tuple.first || (charge == charge_mass_tuple.first && precursor_mass <= charge_mass_tuple.second);
}

bool spectrum::root_scale_intensities() {
    for (float &intensity : intensities) {
        intensity = sqrt(intensity);
    }
    return true;
}

bool spectrum::normalize_intensities() {
    float magnitude = 0.f;
    for (float &i : intensities) {
        magnitude += i*i;
    }
    magnitude = sqrt(magnitude);

    for (float &i : intensities) {
        i = i / magnitude;
    }

    return true;
}

int spectrum::get_mz_bin(float mz) {
    int bin = int((mz - BIN_MIN_MZ) / settings::bin_size);
    return bin;
}

bool spectrum::add_intensity_to_bin(int bin, float intensity) {
    if (bin < 0 || bin > num_bins - 1) {
        return false;
    }
    const vector<int>::iterator &bin_iter = std::find(binned_peaks.begin(), binned_peaks.end(), bin);
    if (bin_iter != binned_peaks.end()) {
        int j = bin_iter - binned_peaks.begin();
        binned_intensities[j] = sqrt(binned_intensities[j] * binned_intensities[j] + intensity * intensity);
    } else {
        binned_peaks.push_back(bin);
        binned_intensities.push_back(intensity);
    }
    return true;
}

