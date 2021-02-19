#include "index_file_writer.h"
#include <fstream>
#include <iostream>

bool index_file_writer::stream_peaks_to_file(std::fstream &f, unsigned int parent_id, const std::shared_ptr<spectrum>& spec) {

    std::string delimiter = ";";
    for (int i = 0; i < spec->binned_peaks.size(); ++i) {
        int bin = spec->binned_peaks[i];
        float intensity = spec->binned_intensities[i];

        f << parent_id << delimiter << bin << delimiter << intensity << "\n";

    }


    return true;
}
