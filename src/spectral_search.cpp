#include <iostream>
#include <cmath>
#include <fstream>
#include "spectral_search.h"
#include "scores.h"

spectral_search::spectral_search() {

}

spectral_search::spectral_search(library *query_lib) : query_lib(query_lib) {

}

spectral_search::spectral_search(library *query_lib, library *target_lib) :
        query_lib(query_lib), target_lib(target_lib)
{

}

bool spectral_search::search_target_library(library *target_lib) {
    this->target_lib = target_lib;
    return search_target_library();
}

bool spectral_search::is_candidate_suitable(spectrum *candidate_spectrum, spectrum *query_spectrum) {
    bool has_equal_charge = candidate_spectrum->charge == query_spectrum->charge;
    bool is_in_mass_range = abs(candidate_spectrum->precursor_mass - query_spectrum->precursor_mass) < mz_tolerance;
    return has_equal_charge && is_in_mass_range;
}

vector<match> spectral_search::get_results() {
    return search_results;
}

bool spectral_search::save_results_to_file(string path, string delimiter) {

    fstream outfile;
    outfile.open(path, ios::out);

    if (!outfile.good())
        return false;

    // Add header
    outfile << "spectrum"+delimiter+"match"+delimiter+"peptide"+delimiter+"dot-product"+delimiter+"mass-difference\n";

    // Go through matches and parse relevant information for each
    for (int i = 0; i < search_results.size(); ++i) {
        match psm = search_results[i];
        outfile << psm.query_spectrum->name << delimiter << psm.matched_spectrum->name << delimiter << psm.matched_spectrum->peptide << delimiter << psm.dot_product << delimiter << psm.mass_difference << endl;
    }

    outfile.close();
    return true;
}

bool spectral_search::read_results_from_file(string path, char delimiter, bool read_dot, bool has_header) {

    /*
     * Requires file format to match delimiter separated result file as it is generated by this program
     * At the very least the format has to respect: name (of query_spectrum) delimiter match (name of library spectrum)
     */

    search_results.clear();
    fstream infile;

    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
    }
    if (read_dot) {
        cerr << "Reading dot-product not implemented. Do so. Proceeding without." << endl;
    }

    string line;
    if (has_header) {
        getline(infile, line);
    }

    string name, match_name;
    while (!infile.eof()) {
        getline(infile, name, delimiter);
        getline(infile, match_name, delimiter);
        if (read_dot) {

        }
        getline(infile, line);

        // Find spectra in the libraries according to the names

        spectrum *query_spectrum, *matched_spectrum;
        for (spectrum *s : query_lib->spectrum_list) { //TODO run-time optimize if this takes to long
            if (s->name == name) {
                query_spectrum = s;
            }
        }
        for (spectrum *s : target_lib->spectrum_list) { //TODO run-time optimize if this takes to long
            if (s->name == match_name) {
                matched_spectrum = s;
            }
        }
        search_results.emplace_back(query_spectrum, matched_spectrum, -1.f, -1);


    }

    return false;
}

bool spectral_search::rescore_matches() {

    for (match &m : search_results) {
        m.dot_product = scores::dot_product(m.query_spectrum->bins, m.matched_spectrum->bins);
    }

    return true;
}

bool spectral_search::search_target_library() {
    cout << "Begin searching target library" << endl;
    search_results.clear();

    for (int i = 0; i < query_lib->spectrum_list.size(); ++i) {
        spectrum* query_spectrum = query_lib->spectrum_list[i];

        if (i % 1000 == 0) {
            cout << "progress: " << i << " of " << query_lib->spectrum_list.size() << " " << (float(i) / query_lib->spectrum_list.size()) * 100 << " %" << endl;
        }

        float max_dot = -1.f;
        spectrum *best_candidate = nullptr;
        //Naive exhaustive spectral_search
        for (spectrum* candidate_spectrum : target_lib->spectrum_list) {
            if (is_candidate_suitable(candidate_spectrum, query_spectrum)) {
                float dot = scores::dot_product(query_spectrum->bins, candidate_spectrum->bins);
                if (dot >= max_dot) {//todo What if equal
                    best_candidate = candidate_spectrum;
                    max_dot = dot;
                }
            }
        }
        if (max_dot >= 0.0) { // if any match was found (i.e. any spectra in mz range)
            match best_match(query_spectrum, best_candidate, max_dot, 1);
            search_results.push_back(best_match);
        }

    }

    return true;
}

