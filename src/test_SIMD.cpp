#include <iostream>
#include <chrono>
#include <vector>
//#include <avxintrin.h>

#include <immintrin.h>
#include <memory>

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
    for (int i = 0; i < vec.size(); i+=8) {
        __m256 _mini_vector = _mm256_load_ps(&vec[i]);
        __m256 _result = _mm256_mul_ps(_scalar, _mini_vector);
        float m[8];
        _mm_store_ps(&m[0], _result);
        _mm256_storeu_ps(&res[i], _result);
    }
}


int main() {

    cout << "Hello SIMD user" << endl;


    /*
     * Args
     */

    int n = 100000000;
    alignas(16) std::vector<float> vec;
    alignas(16) std::vector<float> res;


    cout << alignof(vec) << endl;
    vec.reserve(n);
    res.resize(n);
    //posix_memalign(vec, 16, 16)


    for (int i = 0; i < n; ++i) {
        vec.push_back(0.5f * (i));
    }

    cout << "Test naive:" << endl;
    auto start = chrono::high_resolution_clock::now();
    float val = 0.5f;
    multiply(val, vec, res);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;

    print_res(res);
    res.resize(n, 0);

    cout << "Test m128:" << endl;
    start = chrono::high_resolution_clock::now();
    __m128 _scalar = _mm_set_ps(0.5f, 0.5f, 0.5f, 0.5f);
    multiply128(_scalar, vec, res);

    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;




    print_res(res);
    res.resize(n, 0);

    cout << "Test m256:" << endl;
    start = chrono::high_resolution_clock::now();
    //__m256 _scalar256 = _mm256_set1_ps(0.5f);
    __m256 _scalar256 = _mm256_set_ps(0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f,0.5f);
    multiply256(_scalar256, vec, res);
    stop = chrono::high_resolution_clock::now();
    duration = duration_cast<chrono::microseconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " microseconds" << endl;

    print_res(res);
}
