#include <iostream>
#include <chrono>
#include <vector>

using namespace std;

int main() {

    cout << "Hello SIMD user" << endl;

    auto start = chrono::high_resolution_clock::now();

    /*
     * Args
     */

    std::vector<float> vec;
    vec.reserve(10000);
    for (int i = 0; i < 10000; ++i) {
        vec.push_back(0.5f * i);
    }




    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " seconds" << endl;

}
