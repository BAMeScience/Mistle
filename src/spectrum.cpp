#include <iostream>
#include <cmath>
#include "spectrum.h"
#include "DefineConstants.h"

spectrum::spectrum() {

}

bool spectrum::bin_peaks() {
    int num_bins = BIN_MAX_MZ - BIN_MIN_MZ; // spectrast + 1 // why?

    bins = vector<float>(num_bins, 0.0);

    for (int i = 0; i < peak_positions.size(); ++i) {
        int bin = int(peak_positions[i]) - BIN_MIN_MZ; // this should never be < 0, but should we catch this?
        float intensity = normalize_intensity(intensities[i]); //todo is this may be normalizd already

        bins[bin] = intensity;
    }


    return true;
}

float spectrum::normalize_intensity(float intensity) { //todo cross check spectrast
    return sqrt(intensity);
}
