#include "index_file_writer.h"
#include "DefineConstants.h"
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
    f.precision(FLOAT_OUTPUT_PRECISION);

    //Have num precursors as header (needed for efficient parsing)
    f << "Num: " << precursors.size() << "\n";

    /*
     * ENCODING: ID;RANK;MZ;CHARGE;PEPTIDE
     */
    for (precursor &p : precursors) {
        f << p.id << delimiter << p.rank << delimiter << p.mz << delimiter << p.charge << delimiter << p.peptide << "\n";
    }

    f.close();

    return true;
}

bool index_file_writer::save_matches_to_file(const std::string &file_path, std::vector<match> &matches) {
    std::fstream outfile;
    std::string delimiter = ";";


    outfile.open(file_path, std::ios::out);
    if (!outfile.good())
        return false;

    // Add header
    outfile << "spectrum"+delimiter+"match"+delimiter+"peptide"+delimiter+"dot-product"+delimiter+"mass-difference\n";

    // Go through matches and parse relevant information for each
    for (int i = 0; i < matches.size(); ++i) {
        match psm = matches[i];
        //TODO
        outfile << psm.query_spectrum->name << delimiter << psm.matched_spectrum->name << delimiter << psm.matched_spectrum->peptide << delimiter << psm.dot_product << delimiter << psm.mass_difference << "\n";
    }

    outfile.close();
    return true;
}

bool index_file_writer::stream_peaks_to_binary_file(std::fstream &f, unsigned int parent_id,
                                                    const std::shared_ptr<spectrum> &spec) {
    for (int i = 0; i < spec->intensities.size(); ++i) {
        float mz = spec->peak_positions[i];
        float intensity = spec->intensities[i];

        f.write((char *) &parent_id, sizeof(unsigned int));
        f.write((char *) &mz, sizeof(float));
        f.write((char *) &intensity, sizeof(float));

    }


    return true;
}

bool index_file_writer::save_precursor_index_to_binary_file(const std::string &file_path,
                                                            std::vector<precursor> &precursors) {

    std::fstream f(file_path, std::ios::binary | std::ofstream::out);


    //Save #of precursors as first element (needed for efficient parsing)
    unsigned int size = precursors.size();
    f.write((char *) &size, sizeof(unsigned int));

    /*
     * ENCODING: ID;RANK;MZ;CHARGE;SIZE_OF_PEPTIDE;PEPTIDE (No semi-colons in binary)
     */
    for (precursor &p : precursors) {
        f.write((char *) &p.id, sizeof(unsigned int));
        f.write((char *) &p.rank, sizeof(unsigned int));
        f.write((char *) &p.mz, sizeof(float));
        f.write((char *) &p.charge, sizeof(int));
        size_t pep_size = p.peptide.size();
        f.write((char *) &pep_size, sizeof(pep_size));
        f.write(p.peptide.c_str(), pep_size);
    }

    f.close();

    return true;
}
