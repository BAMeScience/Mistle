#ifndef SIMPLE_EXAMPLE_LIBRARY_H
#define SIMPLE_EXAMPLE_LIBRARY_H
#include <vector>
#include "spectrum.h"

using namespace std;

class library {
public:
    library();
    library(string path);
    bool load_library_from_file(string path);

    vector<spectrum*> spectrum_list;
private:


};


#endif //SIMPLE_EXAMPLE_LIBRARY_H
