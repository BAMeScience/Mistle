

#include <iostream>
#include "indexing_manager.h"

using namespace std;

int main() {
    cout << "Hello World Builder" << endl;


    std::string directory = "/home/ynowatzk/data/9MM/msp";
    indexing_manager im(directory);

    im.build_indices();

}