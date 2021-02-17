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
            lib_files.push_back(entry);
        }
    }


}

bool indexing_manager::build_indices() {
    cout << "Manager too busy to answer phone. Calling back later" << endl;

    for (int i = 0; i < lib_files.size(); ++i) {
        cout << "Parsing library file no. " << i << " (" << lib_files[i].path().filename() << ")" << endl;

        parse_file(lib_files[i].path().string());

    }


    return true;
}

unsigned int indexing_manager::assign_to_index(float mz) {
    for (int i = 0; i < (num_indices - 1); ++i) {
        if (mz < idx_limits[i]) {
            return i;
        }
    }
    return num_indices;
}

bool indexing_manager::parse_file(std::string file_path) {




    return false;
}
