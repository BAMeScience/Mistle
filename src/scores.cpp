#include "scores.h"
#include <cmath>
#include <iostream>

float scores::dot_product(vector<float> &target_bins, vector<float> &other_bins) {
    float dot = 0.f;
    int num_bins = 0;

    for (int i = 0; i < target_bins.size(); ++i) {
        dot += target_bins[i] * other_bins[i];
        /*if (target_bins[i] * other_bins[i] > 0) {
            //cout << target_bins[i] << " * " << other_bins[i] << " = " << target_bins[i] * other_bins[i] << endl;
            ++num_bins;
        }*/
    } //TODO try and compare runtime for iterator

    /*float m1 = 0;
    for (float f:target_bins) {
        m1 += f*f;
    }
    m1 = sqrt(m1);

    float m2 = 0;
    for (float f: other_bins) {
        m2 += f*f;
    }
    m2 = sqrt(m2);

    cout << m1 << " " << m2 << " " << dot << " " << dot / (m1 * m2) << endl;*/
    //cout << " no. " << num_bins << " ";
    return dot;
}
