#include <iostream>
#include <chrono>
#include "scanner.h"

using namespace std;

int main() {

    cout << "Hello Quick Scan" << endl;

    string directory = "/home/ynowatzk/data/9MM/msp/";
    scanner sc;

    cout << "Scanning input:" << endl;
    auto start = chrono::high_resolution_clock::now();
    sc.scan_directory(directory);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Scan Time: " <<  duration.count() << " seconds" << endl;

    sc.analyze();
    sc.print_scan_results();
    sc.save_precursor_distribution_to_file("./precursors.txt");


    return 0;
}