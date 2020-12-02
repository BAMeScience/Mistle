#ifndef SIMPLE_EXAMPLE_SCORING_HANDLER_H
#define SIMPLE_EXAMPLE_SCORING_HANDLER_H
#include <vector>

using namespace std;

class scoring_handler {
    static float dot_product(vector<float> &target_bins, vector<float> &other_bins);
};


#endif //SIMPLE_EXAMPLE_SCORING_HANDLER_H
