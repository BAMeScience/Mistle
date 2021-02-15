#ifndef SIMPLE_EXAMPLE_SCANNER_H
#define SIMPLE_EXAMPLE_SCANNER_H

#include <string>
#include "spectrum.h"

class scanner {
private:

    int no_files_read = 0;
    int kb_lib_size = 0;


public:
    scanner();
    bool scan_directory(std::string path);
    bool scan_file(std::string path);

    bool analyze();
    bool save_precursor_distribution_to_file(std::string path, std::string delimiter="\t");

    bool print_scan_results();

    std::vector<precursor*> parents;
    std::vector<spectrum*> specs;
};


#endif //SIMPLE_EXAMPLE_SCANNER_H
