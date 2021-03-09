#include <iostream>
#include <chrono>
#include "spectrum.h"
#include "msp_reader.h"
#include "scores.h"
#include "spectral_search.h"
#include "library.h"


using namespace std;

int main() {
    cout << "This is a naive search" << endl;
    //FeatureMap fm;
    //Feature feature;
    //PeakSpectrum p;

    //fm.push_back(feature);


    //string msp_file = "/home/ynowatzk/data/pyro_fur/PyroFur_reproduced.msp";
    string msp_file = "/home/ynowatzk/data/9MM/msp/Lactobacillus+acidophilus.msp";
    string mgf_file = "/home/ynowatzk/data/9MM/mgf/9MM_FASP.mgf"; //


    shared_ptr<library> search_lib = make_shared<library>(mgf_file);

    auto start = chrono::high_resolution_clock::now();
    shared_ptr<library> lib = make_shared<library>(msp_file);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;




    /*
     * Search
     */

    start = chrono::high_resolution_clock::now();
    spectral_search search(search_lib, lib);
    search.search_target_library();
    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;

    search.save_results_to_file("./naive_search_results.csv");
    //lib->build_library_index();


    /*
     * SEARCH
     */

    /*
    //Rescoring of spectrast results
    cout << "Reading in" << endl;
    search.read_results_from_file(R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\results\reproduced\sp5_neighbors_matches.tsv)");
    cout << "Rescoring" << endl;
    search.rescore_matches();
    cout << "Saving" << endl;
    search.save_results_to_file(R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\results\reproduced\sp5_neighbors_rescored.tsv)");
    */


    //cout << "Searching fragment ion index" << endl;
    //start = chrono::high_resolution_clock::now();
    //search.search_target_library();
    //stop = chrono::high_resolution_clock::now();
    //duration = duration_cast<chrono::seconds>(stop - start);

    //cout << "Search Time: " <<  duration.count() << " seconds" << endl;
    //search.save_results_to_file("FIIndex.csv");
    //exit(12);

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

