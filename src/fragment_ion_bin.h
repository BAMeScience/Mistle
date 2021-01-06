#ifndef SIMPLE_EXAMPLE_FRAGMENT_ION_BIN_H
#define SIMPLE_EXAMPLE_FRAGMENT_ION_BIN_H
#include <vector>

using namespace std;

struct fragment {
    int parent_id;
    float intensity;
};

class fragment_ion_bin {
public:
    vector<fragment> fragment_list;

    bool score();

};


#endif //SIMPLE_EXAMPLE_FRAGMENT_ION_BIN_H
