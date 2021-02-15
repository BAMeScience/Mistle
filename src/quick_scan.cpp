#include <iostream>
#include <chrono>
#include "scanner.h"
#include "msp_reader.h"
#include "library.h"
#include "spectral_search.h"

using namespace std;

int main() {

    cout << "Hello Quick Scan" << endl;

    string directory = "/home/ynowatzk/data/9MM/msp/Brevibacillus+laterosporus.msp";
    scanner sc;

    cout << "Scanning input:" << endl;
    auto start = chrono::high_resolution_clock::now();
    sc.scan_file(directory);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Scan Time: " <<  duration.count() << " seconds" << endl;

    sc.analyze();
    sc.print_scan_results();
    sc.save_precursor_distribution_to_file("./precursors.txt");


    cout << "Loading spectra from saved positions" << endl;
    start = chrono::high_resolution_clock::now();
    msp_reader::read_spectra_from_positions(directory, sc.parents, sc.specs);
    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;


    // COMPARE SEARCH RESULTS to make sure reading worked
    string mgf_file = "/home/ynowatzk/data/9MM/mgf/9MM_FASP.mgf";
    library *search_lib = new library(mgf_file);
    library *lib = new library(sc.specs);
    lib->build_library_index();

    spectral_search search(search_lib, lib);
    cout << "Searching fragment ion index" << endl;
    start = chrono::high_resolution_clock::now();
    search.search_target_library();
    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::seconds>(stop - start);

    cout << "Search Time: " <<  duration.count() << " seconds" << endl;
    search.save_results_to_file("FIIndex2.csv");


    return 0;
}