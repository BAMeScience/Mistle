#include <iostream>
#include "search_manager.h"

using namespace std;

int main() {
    cout << "Hello World Explorer" << endl;


    std::string search_file = "/home/ynowatzk/data/9MM/mgf/9MM_FASP.mgf";
    std::string index_dir = "./test/";

    search_manager sm(search_file, index_dir);

    sm.prepare_search_library();
    sm.prepare_precursor_index();

    sm.perform_searches();



    /*indexing_manager im(directory);
    auto start = chrono::high_resolution_clock::now();
    im.build_indices();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;
    */



}