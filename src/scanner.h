#ifndef SIMPLE_EXAMPLE_SCANNER_H
#define SIMPLE_EXAMPLE_SCANNER_H

#include <string>
#include "spectrum.h"

using namespace std;

class scanner {
private:

    int no_files_read = 0;
    int kb_lib_size = 0;


public:
    scanner();
    bool scan_directory(string path);
    bool scan_file(string path);

    bool analyze();
    bool save_precursor_distribution_to_file(string path, string delimiter="\t");

    bool print_scan_results();

    vector<precursor*> parents;
    vector<spectrum*> specs;
};


#endif //SIMPLE_EXAMPLE_SCANNER_H
