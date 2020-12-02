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
        //float intensity = rescale_intensity(intensities[i]); //todo how to normalize/rescale??
        float intensity = intensities[i];

        bins[bin] = intensity;
    }


    return true;
}

float spectrum::rescale_intensity(float intensity) { //todo cross check spectrast
    return sqrt(intensity);
}
