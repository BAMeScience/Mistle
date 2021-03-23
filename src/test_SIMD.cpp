#include <iostream>
#include <chrono>
#include <vector>
#include <immintrin.h>

using namespace std;

int main() {

    cout << "Hello SIMD user" << endl;


    /*
     * Args
     */

    int n = 10000000;
    std::vector<float> vec;
    std::vector<float> res;
    vec.reserve(n);
    res.resize(n);
    for (int i = 0; i < n; ++i) {
        vec.push_back(0.5f * (i));
    }


    auto start = chrono::high_resolution_clock::now();
    __m128 _scalar = _mm_set_ps(0.5f, 0.5f, 0.5f, 0.5f);
    for (int i = 0; i < vec.size(); i+=4) {
        __m128 _mini_vector = _mm_loadu_ps(&vec[i]);
        __m128 _result = _mm_mul_ps(_scalar, _mini_vector);
        //float m[4];
        _mm_store_ps(&res[i], _result);
        //cout << m[0] << " " << m[1] << " " << m[2] << " " << m[3] << endl;
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;

    for (int i = 0; i < 10; ++i) {
        cout << res[i] << " ";
    }

    cout << endl;

    for (int i = res.size(); i > (res.size() - 10); --i) {
        cout << res[i-1] << " ";
    }
    cout << endl;

    res.resize(n);
    start = chrono::high_resolution_clock::now();
    float val = 0.5f;
    for (int i = 0; i < vec.size(); i+=1) {
        res[i] = vec[i] * val;
    }

    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;



    for (int i = 0; i < 10; ++i) {
        cout << res[i] << " ";
    }

    cout << endl;

    for (int i = res.size(); i > (res.size() - 10); --i) {
        cout << res[i-1] << " ";
    }
    cout << endl;


}
