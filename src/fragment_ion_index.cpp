#include "fragment_ion_index.h"
#include "DefineConstants.h"

fragment_ion_index::fragment_ion_index(precursor_index *parent_index) {


    fragment_bins = vector<fragment_bin>(BIN_MAX_MZ); //TODO remove/determine actual max #bins
    for (int i = 0; i < parent_index->get_size(); ++i) {
        spectrum *c_spectrum = parent_index->get_spectrum(i);

        /*
         * Iterate all peaks and save them as fragments in the corresponding ion mass bin
         */
        for (int j = 0; j < c_spectrum->binned_peaks.size(); ++j) {
            int bin = c_spectrum->binned_peaks[j];
            fragment frag(c_spectrum->id, c_spectrum->binned_intensities[j]);
            fragment_bins[bin].push_back(frag);
        }
    }

}
