#include "indexing_manager.h"
#include <iostream>
#include <utility>
#include <filesystem>

using namespace std;

indexing_manager::indexing_manager() {

}

indexing_manager::indexing_manager(string path) : path(path) {
    cout << "Calling Manager" << endl;

    for (const auto & entry : std::filesystem::directory_iterator(path)) {
        if (entry.path().extension() == ".msp") {
            cout << entry.path().string() << endl;
            lib_files.push_back(entry);
        }
        //cout << spectrum_list.size() << endl;
    }


}

bool indexing_manager::build_indices() {
    cout << "Manager too busy to answer phone. Calling back later" << endl;


    return true;
}
