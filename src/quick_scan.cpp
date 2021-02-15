#include <iostream>
#include <chrono>
#include "scanner.h"
#include "msp_reader.h"

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

    return 0;
}