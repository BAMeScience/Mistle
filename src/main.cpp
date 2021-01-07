#include <iostream>
#include <numeric>
#include <cmath>
#include <chrono>
#include "spectrum.h"
#include "msp_reader.h"
#include "scores.h"
#include "spectral_search.h"
#include "library.h"


int main() {
    cout << "Welcome, welcome" << endl;


    string msp_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_Complete_simulatedSpectra\PyroFur_reproduced.msp)";
    string mgf_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_SearchFile\pfu_velos.mgf)";
    //vector<spectrum*> library = msp_reader::read_file(R"(C:\Users\ynowatzk\Desktop\data\9MM\simulated_spectra\Brevibacillus+laterosporus.msp)");


    library *search_lib = new library(mgf_file);
    library *lib = new library(msp_file);
    lib->build_library_index();


    /*
     * SEARCH
     */

    spectral_search search(search_lib, lib);

    /*
    //Rescoring of spectrast results
    cout << "Reading in" << endl;
    search.read_results_from_file(R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\results\reproduced\sp5_neighbors_matches.tsv)");
    cout << "Rescoring" << endl;
    search.rescore_matches();
    cout << "Saving" << endl;
    search.save_results_to_file(R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\results\reproduced\sp5_neighbors_rescored.tsv)");
    */


    cout << "Searching fragment ion index" << endl;
    auto start = chrono::high_resolution_clock::now();
    search.search_target_library();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);

    cout << "Search Time: " <<  duration.count() << " seconds" << endl;
    search.save_results_to_file("FIIndex.csv");
    exit(12);

    /*
    auto start = chrono::high_resolution_clock::now();
    search.search_target_library();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);

    cout << "Search Time: " <<  duration.count() << " seconds" << endl;


    vector<match> matches = search.get_results();

    for (int i = 0; i < 10; ++i) {
        cout << matches[i].query_spectrum->name << " " << matches[i].matched_spectrum->peptide << " " << matches[i].dot_product << endl;
    }

    search.save_results_to_file("./out.csv");

    */
    return 0;
}
