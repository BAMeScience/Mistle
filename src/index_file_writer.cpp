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

bool index_file_writer::save_precursor_index(const std::string& file_path, std::vector<precursor> &precursors) {

    std::ofstream f(file_path, std::ofstream::out);
    std::string delimiter = ";";

    //Have num precursors as header (needed for efficient parsing)
    f << "Num: " << precursors.size();

    /*
     * ENCODING: ID;RANK;MZ;CHARGE;PEPTIDE
     */
    for (precursor &p : precursors) {
        f << p.id << delimiter << p.rank << delimiter << p.mass << delimiter << p.charge <<  delimiter << p.peptide << "\n";
    }

    f.close();

    return true;
}
