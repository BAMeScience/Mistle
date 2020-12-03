#include "scores.h"

float scores::dot_product(vector<float> &target_bins, vector<float> &other_bins) {
    float dot = 0;

    for (int i = 0; i < target_bins.size(); ++i) {
        dot += target_bins[i] * other_bins[i];
    } //TODO try and compare runtime for iterator

    return dot;
}
