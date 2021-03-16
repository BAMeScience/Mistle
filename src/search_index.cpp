#include <iostream>
#include <chrono>
#include <utility>
#include "search_manager.h"
#include "thread_pool.h"

using namespace std;
/*
template<typename F, typename... Args>
void time_function(F func, Args&&... args){
    auto start = chrono::high_resolution_clock::now();
    func(std::forward<Args>(args)...);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Time elapsed: " <<  duration.count() << " seconds" << endl;
}

auto time_func =
        [](auto&& func, auto&&... params) {
            // get time before function invocation
            auto start = chrono::high_resolution_clock::now();
            std::forward<decltype(func)>(func)(std::forward<decltype(params)>(params)...);
            auto stop = chrono::high_resolution_clock::now();
            auto duration = duration_cast<chrono::seconds>(stop - start);
            cout << "Time elapsed: " <<  duration.count() << " seconds" << endl;
        };
*/

int main() {

    cout << "Hello World Explorer" << endl;

    thread_pool pool(4);
    std::vector<int> test;

    pool.enqueue([&](){
        pool.mtx.lock();
        test.push_back(1);
        test.push_back(2);
        test.push_back(3);
        test.push_back(4);
        test.push_back(5);
        test.push_back(6);
        test.push_back(1);
        test.push_back(2);
        test.push_back(3);
        test.push_back(4);
        test.push_back(5);
        test.push_back(6);
        pool.mtx.unlock();
        cout << "Hello Task " << endl;
    });
    pool.enqueue([&]() {
        pool.mtx.lock();
        test.push_back(1);
        test.push_back(2);
        test.push_back(3);
        test.push_back(4);
        test.push_back(5);
        test.push_back(6);
        test.push_back(1);
        test.push_back(2);
        test.push_back(3);
        test.push_back(4);
        test.push_back(5);
        test.push_back(6);
        pool.mtx.unlock();
        cout << "Hello Again " << endl;
    });




    auto start = chrono::high_resolution_clock::now();

    /*
     * Args
     */

    std::string search_file = "/home/ynowatzk/data/9MM/mgf/9MM_FASP.mgf";
    std::string index_dir = "./test/";


    /*
     * Preparation and Search
     */

    search_manager sm(search_file, index_dir);

    std::cout << "Preparing libraries and indices" << std::endl;
    sm.prepare_search_library();
    sm.prepare_precursor_index();

    for (int &i : test) {
        cout << " " << i;
    }
    cout << endl;

    exit(12);

    std::cout << "Searching fragment-ion-indices in batches" << endl;
    sm.perform_searches();
    std::cout << "Merging overlapping results" << std::endl;
    auto check_point = chrono::high_resolution_clock::now();
    sm.merge_matches();
    std::cout << "Time elapsed: " << duration_cast<chrono::seconds>(check_point - chrono::high_resolution_clock::now()).count() << " seconds" << endl;
    std::cout << "Writing results to file" << std::endl;
    sm.save_search_results_to_file(index_dir + "results.csv");




    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " seconds" << endl;


    /*indexing_manager im(directory);
    auto start = chrono::high_resolution_clock::now();
    im.build_indices();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;
    */



}