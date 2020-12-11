#include <iostream>
#include <cmath>
#include <numeric>
#include "spectrum.h"
#include "DefineConstants.h"

spectrum::spectrum() {

}

bool spectrum::bin_peaks(bool root_rescale, bool normalize) {
    int num_bins = BIN_MAX_MZ - BIN_MIN_MZ; // spectrast + 1 // why?

    bins = vector<float>(num_bins, 0.0);

    for (int i = 0; i < peak_positions.size(); ++i) {
        int bin = int(peak_positions[i]) - BIN_MIN_MZ;
        if (bin < 0 || bin > BIN_MAX_MZ) { // TODO spectraST light ion, cut_off. bin < 180 ??
            //TODO what would spectrast do?
            //cout << "Warning peak out of bin range :: discarding intensity" << endl;
            //bin = 0;

            continue;
        }
        if (abs(float(bin) - precursor_mass) < 10)
            continue;

        float intensity = intensities[i];
        if (root_rescale)
            intensity = sqrt(intensity);

        bins[bin] += intensity; //+= accumulate if multiple peaks fall into the sam bin
        if (intensity_bin_spanning_factor > 0.f) {
            float neighbor_intensity = intensity * intensity_bin_spanning_factor;
            if (bin > 0) {
                bins[bin-1] += neighbor_intensity;
            }
            if (bin < BIN_MAX_MZ) {
                bins[bin+1] += neighbor_intensity;
            }
        }
    }

    if (normalize) {
        /*normalize_bins();
        float min_cut = 0.05f; //TODO undirty or delete this
        for (float &i : bins) {
            if (i<min_cut) {
                i = 0.f;
            }
        }*/
        return normalize_bins();
    }
    return true;
}

float spectrum::rescale_intensity(float intensity) { //todo cross check spectraST
    return sqrt(intensity);
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
