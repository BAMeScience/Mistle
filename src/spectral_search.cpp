#include <iostream>
#include <cmath>
#include <fstream>
#include "spectral_search.h"
#include "scores.h"

spectral_search::spectral_search() {

}

spectral_search::spectral_search(library *search_lib) : search_lib(search_lib) {

}

bool spectral_search::search_target_library(library *target_lib) {
    this->target_lib = target_lib;
    cout << "Begin searching target library" << endl;


    for (int i = 0; i < search_lib->spectrum_list.size(); ++i) {
        spectrum* query_spectrum = search_lib->spectrum_list[i];

        if (i % 1000 == 0) {
            cout << "progress: " << i << " of " << search_lib->spectrum_list.size() << " " << (float(i) / search_lib->spectrum_list.size()) * 100 << " %" << endl;
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

vector<match> spectral_search::search_specific_by_name(string query_name, library *target_lib, int num_matches) {

    spectrum *query_spectrum = nullptr;
    for (spectrum *s : search_lib->spectrum_list) {
        if (s->name == query_name) {
            query_spectrum = s;
            break;
        }
    }
    if (query_spectrum == nullptr) {
        cerr << "no spectrum find under such name" << endl;
        return vector<match>();
    }

    //TODO w.i.p What am I even doing here?
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




    return vector<match>();
}

