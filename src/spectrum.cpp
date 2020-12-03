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
    float magnitude = 0;

    for (int i = 0; i < peak_positions.size(); ++i) {
        int bin = int(peak_positions[i]) - BIN_MIN_MZ; // this should never be < 0, but should we catch this?
        if (bin < BIN_MIN_MZ || bin > BIN_MAX_MZ) {
            //TODO what would spectrast do?
            //cout << "Warning peak out of bin range :: discarding intensity" << endl;
            continue;
        }
        //float intensity = rescale_intensity(intensities[i]); //todo how to normalize/rescale??
        float intensity = intensities[i];
        if (root_rescale)
            intensity = sqrt(intensity);

        bins[bin] += intensity; //+= accumulate if multiple peaks fall into the sam bin
        magnitude += (intensity * intensity);
    }

    if (normalize) {
        return normalize_bins(sqrt(magnitude));
    }
    return true;
}

float spectrum::rescale_intensity(float intensity) { //todo cross check spectraST
    return sqrt(intensity);
}

bool spectrum::normalize_bins(float magnitude) {
    if (magnitude < 0) {
        magnitude = 0;
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
