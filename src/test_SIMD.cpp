#include <iostream>
#include <chrono>
#include <vector>
//#include <avxintrin.h>

#include <immintrin.h>
#include <memory>
#include <complex>
#include "library.h"

using namespace std;

void print_res(std::vector<float> &res) {
    for (int i = 0; i < 10; ++i) {
        cout << res[i] << " ";
    }
    cout << endl;
    for (int i = res.size(); i > (res.size() - 10); --i) {
        cout << res[i-1] << " ";
    }
    cout << endl;
}

void multiply(float scalar, std::vector<float> &vec, std::vector<float> &res) {
    for (int i = 0; i < vec.size(); ++i) {
        res[i] = scalar * vec[i];
    }
}

void multiply128(__m128 _scalar, std::vector<float> &vec, std::vector<float> &res) {
    for (int i = 0; i < vec.size(); i+=4) {
        __m128 _mini_vector = _mm_load_ps(&vec[i]);
        __m128 _result = _mm_mul_ps(_scalar, _mini_vector);
        //float m[4];
        _mm_store_ps(&res[i], _result);
    }

}

void multiply256(__m256 _scalar, std::vector<float> &vec, std::vector<float> &res) {
    for (int i = 0; i < vec.size() - 100; i+=8) {
        __m256 _mini_vector = _mm256_load_ps(&vec[i + 2]);
        __m256 _result = _mm256_mul_ps(_scalar, _mini_vector);

        _mm256_storeu_ps(&res[i], _result);
    }
}

float normal_pdf(float x, float m, float s)
{
    static const float inv_sqrt_2pi = 0.3989422804014327;
    float a = (x - m) / s;

    return inv_sqrt_2pi / s * std::exp(-0.5f * a * a);
}

float normal_pdf_div(float x, float m, float s)
{
    return normal_pdf(x,m,s) / normal_pdf(m,m,s);
}

int main() {

    library lib;
    lib.is_indexed;

    cout << "Hello SIMD user" << endl;

    for (int i = -3; i <= 3; ++i) {
        cout << i << ":\t" << normal_pdf_div(i,0,1) << endl;
    }
    float i1 = 0.3f;
    float i2 = 0.25f;
    float dist = 0.1; //Daltons
    float s = 0.1;

    cout << i1 * i2  << " " << i1 * normal_pdf_div(dist, 0,s) * i2 << " " << i1 * normal_pdf_div(dist/2.f, 0, s)  * normal_pdf_div(dist/2, 0, s) * i2 << endl;
    exit(12);
    /*
     * Args
     */
    volatile int8_t test;
    cout << test << endl;
    int n = 1000;
    alignas(32) std::vector<std::vector<float>> vec(2);
    alignas(32) std::vector<float> res;


    cout << alignof(vec) << endl;
    cout << alignof(vec[1]) << endl;
    vec[1].reserve(n);
    res.resize(n);
    //posix_memalign(vec, 16, 16)


    for (int i = 0; i < n; ++i) {
        vec[1].push_back(0.5f * (i));
    }

    cout << "Test naive:" << endl;
    auto start = chrono::high_resolution_clock::now();
    float val = 0.5f;
    multiply(val, vec[1], res);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;

    print_res(res);
    res.resize(n, 0);

    cout << "Test m128:" << endl;
    start = chrono::high_resolution_clock::now();
    __m128 _scalar = _mm_set_ps(0.5f, 0.5f, 0.5f, 0.5f);
    multiply128(_scalar, vec[1], res);

    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;




    print_res(res);
    res.resize(n, 0);

    cout << "Test m256:" << endl;
    start = chrono::high_resolution_clock::now();
    //__m256 _scalar256 = _mm256_set1_ps(0.5f);
    __m256 _scalar256 = _mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f);
    multiply256(_scalar256, vec[1], res);
    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;

    print_res(res);
}
